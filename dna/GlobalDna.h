#ifndef GLOBAL_DNA_H
#define GLOBAL_DNA_H

#include <vector>
#include <iostream>
using namespace std;

#define u8 unsigned char
#define u32 unsigned int
#define u64 unsigned long long


//================================= VINT BEGIN =================================
u8 vint_buf[5]; //buffer for holding the vlength-bytes of a number
u8 vint_bytes; //number of bytes

void append_vint(vector<u8> & collector) //set number_buf & number_bytes properly first
{
	for(int i=0; i<vint_bytes; i++) collector.push_back(vint_buf[i]);
}

//-----------------

void to_vint(u32 i)
{
	if(i == 0)
	{
		vint_buf[0] = 0x80;
		vint_bytes = 1;
	}
	else //i > 0
	{
		vint_bytes = 0;
		while(i)
		{
			vint_buf[vint_bytes] = (i & 0x7f);
			vint_bytes ++;
			i >>= 7;
		}
		vint_buf[vint_bytes - 1] |= 0x80;
	}
}

u32 parse_vint(u8 * & head)
{
	int i=0;
	u32  ret = 0;
	while(!(head[i] & 0x80))
	{
		ret += (head[i] << 7*i);
		i++;
	}
	ret += ((head[i] & 0x7f) << 7*i);
	head += i + 1;
	return ret;
}

void parse_vints(vector<u32> & collector, vector<u8> & vec_vints)
{
	u8 * head = &vec_vints[0];
	u8 * p = head;
	u32  size = vec_vints.size();
	while(p - head < size)
	{
		u32 cur = parse_vint(p);
		collector.push_back(cur);
	}
}
//================================= VINT END =================================


//================================= K-MER BEGIN =================================
typedef u64 int k_mer;

//mapping:
// A: 00
// T: 11
// G: 10
// C: 01
static const k_mer A = 0ull;
static const k_mer T = 3ull;
static const k_mer G = 2ull;
static const k_mer C = 1ull;

// 0, 1, 2, 3 (A, C, G, T) -> 0x08, 0x01, 0x02, 0x04
static const u32 bit_pos[4] = {8u, 1u, 2u, 4u};

//====================================
//bits to edges mapping (in one byte):
static const u8 inA = 0x80;
static const u8 inT = 0x40;
static const u8 inG = 0x20;
static const u8 inC = 0x10;
static const u8 outA = 0x08;
static const u8 outT = 0x04;
static const u8 outG = 0x02;
static const u8 outC = 0x01;

static const u8 ATGC_bits[8] = {inA, inT, inG, inC, outA, outT, outG, outC};
static const k_mer bit2ATGC[8] = {A, T, G, C, A, T, G, C};

int mer_length;
k_mer kick;
k_mer END_MER;
k_mer NULL_MER = 0x8000000000000000; //k-mer does not use the highest 2 bits

void set_mer_length(int len)
{
	mer_length = len;
	kick = 0xFFFFFFFFFFFFFFFF >> (64 - 2 * mer_length);
}

void convert(k_mer val, char* buf, int length)
{
	for(int i = 0; i < length; i++)
	{
		k_mer tag = (val & 3ull);
		if(tag == T) buf[length-1-i] = 'T';
		else if(tag == G) buf[length-1-i] = 'G';
		else if(tag == C) buf[length-1-i] = 'C';
		else buf[length-1-i] = 'A';
		val >>= 2;
	}
}

k_mer getRC(k_mer id)
{
	k_mer rc = 0;
	for(int i = 0; i < mer_length; i++)
	{
		k_mer tag =  id & 3ull; //last 2 bits
		tag ^= 3ull; //complement
		rc |= tag; //append to rc
		id >>= 2; //next 2 bits
		if(i != mer_length - 1) rc <<= 2;
	}
	return rc;
}

int getShift(bool v1_pol, bool v2_pol)//four bytes, each for LL, LH, HL, and HH
{
	int shift = 0;
	if (!v1_pol) shift += 16;
	if (!v2_pol) shift += 8;
	return shift;
}

k_mer get_neighbor(k_mer id, u8 bit, bool v1_pol, bool v2_pol)  //id is the ID of the current vertex
{
	k_mer neighborID = id;
	//------
	if (bit == inA)
	{
		if(v2_pol) neighborID = getRC(id);
		neighborID >>= 2;
		if(v1_pol) neighborID = getRC(neighborID);
	}
	else if (bit == inT)
	{
		if(v2_pol) neighborID = getRC(id);
		neighborID >>= 2;
		neighborID |= (T << (2 * mer_length - 2));
		if(v1_pol) neighborID = getRC(neighborID);
	}
	else if (bit == inG)
	{
		if(v2_pol) neighborID = getRC(id);
		neighborID >>= 2;
		neighborID |= (G << (2 * mer_length - 2));
		if(v1_pol) neighborID = getRC(neighborID);
	}
	else if (bit == inC)
	{
		if(v2_pol) neighborID = getRC(id);
		neighborID >>= 2;
		neighborID |= (C << (2 * mer_length - 2));
		if(v1_pol) neighborID = getRC(neighborID);
	}
	else if (bit == outA)
	{
		if(v1_pol) neighborID = getRC(id);
		neighborID <<= 2;
		neighborID &= kick; //kick out highest two bits
		if(v2_pol) neighborID = getRC(neighborID);
	}
	else if (bit == outT)
	{
		if(v1_pol) neighborID = getRC(id);
		neighborID <<= 2;
		neighborID &= kick;
		neighborID |= T;
		if(v2_pol) neighborID = getRC(neighborID);
	}
	else if (bit == outG)
	{
		if(v1_pol) neighborID = getRC(id);
		neighborID <<= 2;
		neighborID &= kick;
		neighborID |= G;
		if(v2_pol) neighborID = getRC(neighborID);
	}
	else if (bit == outC)
	{
		if(v1_pol) neighborID = getRC(id);
		neighborID <<= 2;
		neighborID &= kick;
		neighborID |= C;
		if(v2_pol) neighborID = getRC(neighborID);
	}
	return neighborID;
}

//special bitmap for dead end
u8 DEAD_END = 0x80;

struct neighbor_info
{
	u8 bitmap;
	//- - - A A B C C
	//AA: for ATGC
	//B: for in/out, in = 1, out = 0
	//CC: for L/H:L/H, L = 0, H = 1

	friend ibinstream& operator<<(ibinstream& m, const neighbor_info& v)
	{
		m << v.bitmap;
		return m;
	}

	friend obinstream& operator>>(obinstream& m, neighbor_info& v)
	{
		m >> v.bitmap;
		return m;
	}

	//------

	neighbor_info()
	{
		bitmap = 0;
	}

	k_mer get_ATGC()
	{
		return (bitmap >> 3);
	}

	bool is_in()
	{
		return (bitmap & 0x04);
	}

	bool left_isH()
	{
		return (bitmap & 0x02);
	}

	bool right_isH()
	{
		return (bitmap & 0x01);
	}

	//------

	void set_ATGC(k_mer tag)
	{
		bitmap &= 0xE7;//eraser mask: 11100111
		bitmap |= (tag << 3);
	}

	void set_in(bool tag)
	{
		if(tag) bitmap |= 0x04;
		else bitmap &= 0xFB;//eraser mask: 11111011
	}

	void set_left(bool isH)
	{
		if(isH) bitmap |= 0x02;
		else bitmap &= 0xFD;//eraser mask: 11111101
	}

	void set_right(bool isH)
	{
		if(isH) bitmap |= 0x01;
		else bitmap &= 0xFE;//eraser mask: 11111110
	}

	//------

	bool dead_end()
	{
		return bitmap == DEAD_END;
	}

	void reverse()
	{
		bitmap ^= 0x1F; //reverse last 5 bits
		//swap last 2 bits
		u8 low = bitmap & 0x01;
		u8 high = bitmap & 0x02;
		bitmap &= 0xFC; //clear last two bits
		bitmap |= (low << 1);
		bitmap |= (high >> 1);
	}
};

k_mer get_neighbor(k_mer me, neighbor_info nb)
{
	k_mer tag = nb.get_ATGC();
	tag = bit_pos[tag];
	if(nb.is_in()) tag <<= 4;
	return get_neighbor(me, tag, nb.left_isH(), nb.right_isH());
}

//========================== k-mer pair ==========================
struct kmer_pair
{
	k_mer v1;
	k_mer v2;

	kmer_pair()
	{
	}

	kmer_pair(k_mer v1, k_mer v2)
	{
		if(v1 < v2)
		{
			this->v1 = v1;
			this->v2 = v2;
		}
		else
		{
			this->v1 = v2;
			this->v2 = v1;
		}
	}

	bool set(k_mer v1, k_mer v2)
	{
		if(v1 < v2)
		{
			this->v1 = v1;
			this->v2 = v2;
			return true;
		}
		else
		{
			this->v1 = v2;
			this->v2 = v1;
			return false;
		}
	}

	inline bool operator<(const kmer_pair& rhs) const
	{
		return (v1 < rhs.v1) || ((v1 == rhs.v1) && (v2 < rhs.v2));
	}

	inline bool operator>(const kmer_pair& rhs) const
	{
		return (v1 > rhs.v1) || ((v1 == rhs.v1) && (v2 > rhs.v2));
	}

	inline bool operator==(const kmer_pair& rhs) const
	{
		return (v1 == rhs.v1) && (v2 == rhs.v2);
	}

	inline bool operator!=(const kmer_pair& rhs) const
	{
		return (v1 != rhs.v1) || (v2 != rhs.v2);
	}

	int hash()
	{
		size_t seed = 0;
		hash_combine(seed, v1);
		hash_combine(seed, (v1 >> 32));
		hash_combine(seed, v2);
		hash_combine(seed, (v2 >> 32));
		return seed % ((unsigned int)_num_workers);
	}
};

ibinstream& operator<<(ibinstream& m, const kmer_pair& v)
{
	m << v.v1;
	m << v.v2;
	return m;
}

obinstream& operator>>(obinstream& m, kmer_pair& v)
{
	m >> v.v1;
	m >> v.v2;
	return m;
}

class KmerPairHash
{
public:
	inline int operator()(kmer_pair& key)
	{
		return key.hash();
	}
};

namespace __gnu_cxx
{
template <>
struct hash<kmer_pair>
{
	size_t operator()(kmer_pair& pair) const
	{
		size_t seed = 0;
		hash_combine(seed, pair.v1);
		hash_combine(seed, (pair.v1 >> 32));
		hash_combine(seed, pair.v2);
		hash_combine(seed, (pair.v2 >> 32));
		return seed;
	}
};
}

//=======================================

struct ATGC_bitmap
{
	vector<u8> seq;
	int length;

	ATGC_bitmap()
	{
		length = 0;
	}

	void init(k_mer mer)
	{
		length = mer_length;
		//------
		int byte_num = mer_length / 4;
		int mod = mer_length % 4;
		if(mod) byte_num++;
		mer <<= (64 - 2*mer_length);//align to right
		for(int i=0; i<byte_num; i++)
		{
			k_mer byte = mer >> 8*(7-i);
			byte &= 0xFF;
			seq.push_back(byte);
		}
	}

	void append(k_mer tag)
	{
		length += 1;
		//decide whether need to increase "seq"
		int byte_num = length / 4;
		int mod = length % 4;
		if(mod) byte_num++;
		if(byte_num > seq.size() || seq.size() == 0) seq.push_back(0);
		//set the bits
		if(mod == 0) seq.back() |= (u8)tag; //00 00 00 11
		else if(mod == 1) seq.back() |= (u8)(tag << 6); //11 00 00 00
		else if(mod == 2) seq.back() |= (u8)(tag << 4); //00 11 00 00
		else seq.back() |= (u8)(tag << 2); //00 00 11 00
	}

	k_mer get(int pos)
	{
		int byte_no = pos / 4;
		int mod = pos % 4;
		if(mod == 0) return (seq[byte_no] & 0xC0) >> 6;//11 00 00 00
		else if(mod == 1) return (seq[byte_no] & 0x30) >> 4;//00 11 00 00
		else if(mod == 2) return (seq[byte_no] & 0x0C) >> 2;//00 00 11 00
		else return (seq[byte_no] & 0x03);//00 00 00 11
	}

	ATGC_bitmap get_reverse()
	{
		ATGC_bitmap rev;
		for(int i = 0; i < length; i++)
		{
			rev.append(get(length - i -1) ^ 3ull);
		}
		rev.length = length;
		return rev;
	}

	string toString()
	{
		string re;
		for(int i=0; i<length; i++)
		{
			k_mer cur = get(i);
			if(cur == A) re += "A";
			else if(cur == T) re += "T";
			else if(cur == G) re += "G";
			else if(cur == C) re += "C";
		}
		return re;
	}

	friend ibinstream& operator<<(ibinstream& m, const ATGC_bitmap& v)
	{
		m << v.seq;
		m << v.length;
		return m;
	}

	friend obinstream& operator>>(obinstream& m, ATGC_bitmap& v)
	{
		m >> v.seq;
		m >> v.length;
		return m;
	}
};

int edist(ATGC_bitmap bitmap1, ATGC_bitmap bitmap2)
{
	int rows = bitmap1.length + 1;
	int cols = bitmap2.length + 1;
	int ** d = new int * [rows];
	//------
	for(int i=0; i<rows; i++)
	{
		d[i] = new int[cols];
		d[i][0] = i;
	}
	for(int j=0; j<cols; j++) d[0][j] = j;
	//------
	for(int i=1; i<rows; i++)
	{
		u8 ci = bitmap1.get(i-1);
		for(int j=1; j<cols; j++)
		{
			u8 cj = bitmap2.get(j-1);
			if(ci == cj) d[i][j] = d[i-1][j-1];
			else
			{
				int min = d[i-1][j-1];
				if(min > d[i][j-1]) min = d[i][j-1];
				if(min > d[i-1][j]) min = d[i-1][j];
				d[i][j] = min + 1;
			}
		}
	}
	int result = d[rows-1][cols-1];
	//------
	for(int i=0; i<rows; i++) delete[] d[i];
	delete[] d;
	return result;
}

//=======================================

class Contig
{
public:
	k_mer id; //set top bit to one, different from normal vertices
	ATGC_bitmap seq;
	k_mer in_neighbor; //may be NULL_MER, then in_pol is invalid
	bool in_pol;//(amb) -> [?]:L -> (contig)
	u32 in_count;
	k_mer out_neighbor; //may be NULL_MER, then out_pol is invalid
	bool out_pol;//(contig) -> L:[?] -> (amb)
	u32 out_count;
	int freq;
	bool mark;

	void parse(char* line)
	{
		char * pch;
		pch=strtok(line, "\t");
		id = strtoull(pch, NULL, 10);
		pch=strtok(NULL, " ");
		in_neighbor= strtoull(pch, NULL, 10);
		pch=strtok(NULL, " ");
		in_pol = atoi(pch);
		pch=strtok(NULL, " ");
		in_count = atoi(pch);
		pch=strtok(NULL, " ");
		out_neighbor = strtoull(pch, NULL, 10);
		pch=strtok(NULL, " ");
		out_pol = atoi(pch);
		pch=strtok(NULL, " ");
		out_count = atoi(pch);
		pch=strtok(NULL, " ");
		freq = atoi(pch);
		pch=strtok(NULL, " ");
		seq.length = atoi(pch);
		pch=strtok(NULL, " ");
		int size = atoi(pch);
		for(int i = 0 ; i < size; i++)
		{
			pch=strtok(NULL, " ");
			seq.seq.push_back((u8)atoi(pch));
		}
	}

	void dumpTo(BufferedWriter* writer)
	{
		char buf[100];
		writer->check();
		sprintf(buf, "%llu\t", id);
		writer->write(buf);
		//------
		if(in_neighbor == NULL_MER) sprintf(buf, "%llu 0 0", in_neighbor);
		else sprintf(buf, "%llu %d %u", in_neighbor, (in_pol ? 1 : 0), in_count);
		writer->write(buf);
		//------
		if(out_neighbor == NULL_MER) sprintf(buf, " %llu 0 0", out_neighbor);
		else sprintf(buf, " %llu %d %u", out_neighbor, (out_pol ? 1 : 0), out_count);
		writer->write(buf);
		//------
		sprintf(buf, " %u %u %u", freq, seq.length, seq.seq.size());
		writer->write(buf);
		//------
		for(int i=0; i<seq.seq.size(); i++)
		{
			sprintf(buf, " %u", seq.seq[i]);
			writer->write(buf);
		}
		writer->write("\n");
	}

	friend ibinstream& operator<<(ibinstream& m, const Contig& v)
	{
		m << v.id;
		m << v.seq;
		m << v.in_neighbor;
		m << v.in_pol;
		m << v.in_count;
		m << v.out_neighbor;
		m << v.out_pol;
		m << v.out_count;
		m << v.freq;
		m << v.mark;
		return m;
	}

	friend obinstream& operator>>(obinstream& m, Contig& v)
	{
		m >> v.id;
		m >> v.seq;
		m >> v.in_neighbor;
		m >> v.in_pol;
		m >> v.in_count;
		m >> v.out_neighbor;
		m >> v.out_pol;
		m >> v.out_count;
		m >> v.freq;
		m >> v.mark;
		return m;
	}
};

//=======================================
struct AmbiNB
{
	neighbor_info nb_info;
	int count;
};

ibinstream & operator<<(ibinstream & m, const AmbiNB & v)
{
	m<<v.nb_info;
	m<<v.count;
	return m;
}

obinstream & operator>>(obinstream & m, AmbiNB & v)
{
	m>>v.nb_info;
	m>>v.count;
	return m;
}

struct ContigNB
{
	k_mer nid; //ID of the neighbor on the other end of a contig
	neighbor_info ninfo; //only use in/out, L/H, L/H
	k_mer contigID; //ID of the contig
	int length; //length of the contig
};

ibinstream & operator<<(ibinstream & m, const ContigNB & v)
{
	m << v.nid;
	m << v.ninfo;
	m << v.contigID;
	m << v.length;
	return m;
}

obinstream & operator>>(obinstream & m, ContigNB & v)
{
	m >> v.nid;
	m >> v.ninfo;
	m >> v.contigID;
	m >> v.length;
	return m;
}

//v-type:
static const u8 V_1 =  1;
static const u8 V1_1 = 2;
static const u8 Vm_n = 3;

#endif
