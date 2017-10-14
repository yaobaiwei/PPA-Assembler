#ifndef DNAMESSAGEBUFFER_H_
#define DNAMESSAGEBUFFER_H_

#include "utils/global.h"
#include "utils/communication.h"
#include "utils/vecs.h"
using namespace std;

template <class VertexT, class KeyType, class MessageType, class HashType>
class DNAMessageBuffer
{
public:
	typedef  KeyType KeyT;
	typedef  MessageType MessageT;
	typedef  HashType HashT;
	typedef vector<MessageT> MessageContainerT;
	typedef hash_map<KeyT, int> Map; //int = position in v_msg_bufs //CHANGED FOR VADD
	typedef Vecs<KeyT, MessageT, HashT> VecsT;
	typedef typename VecsT::Vec Vec;
	typedef typename VecsT::VecGroup VecGroup;
	typedef typename Map::iterator MapIter;

	VecsT out_messages;
	Map in_messages;
	vector<MessageContainerT> v_msg_bufs;

	void init(vector<VertexT*> & vertexes)
	{
		v_msg_bufs.resize(vertexes.size());
		for (int i = 0; i < vertexes.size(); i++)
		{
			VertexT* v = vertexes[i];
			in_messages[v->id] = i; //CHANGED FOR VADD
		}
	}

	void add_message(const KeyT& id, const MessageT& msg)
	{
		out_messages.append(id, msg);
	}

	void sync_messages()
	{
		int np = get_num_workers();
		//exchange msgs
		all_to_all(out_messages.getBufs());

		//================================================
		// gather all messages
		for (int i = 0; i < np; i++)
		{
			Vec& msgBuf = out_messages.getBuf(i);
			for (int j = 0; j < msgBuf.size(); j++)
			{
				MapIter it = in_messages.find(msgBuf[j].key);
				if (it != in_messages.end()) //filter out msgs to non-existent vertices
					v_msg_bufs[it->second].push_back(msgBuf[j].msg); //CHANGED FOR VADD
			}
		}
		out_messages.clear();
	}

	vector<MessageContainerT>& get_v_msg_bufs()
	{
		return v_msg_bufs;
	}
};

#endif /* DNAMESSAGEBUFFER_H_ */
