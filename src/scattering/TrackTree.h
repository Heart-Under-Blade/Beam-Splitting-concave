#pragma once

#include <vector>

struct TrackNode
{
	bool isCompleted;
	int facetId;
	TrackNode* m_parent;
	std::vector<TrackNode*> children;

	TrackNode(int value, TrackNode* parent = nullptr)
	{
		isCompleted = false;
		facetId = value;
		m_parent = parent;
	}

	~TrackNode()
	{
		for (auto child : children)
		{
			delete child;
		}
	}

	TrackNode* AddChild(int value)
	{
		TrackNode* child = new TrackNode(value, this);
		children.push_back(child);
		return child;
	}

	TrackNode *FindNode(int value)
	{
		for (auto child : children)
		{
			if (child->facetId == value)
			{
				return child;
			}
		}

		return nullptr;
	}
};

class TrackTree
{
public:
	TrackTree();

private:
	TrackNode *m_root;
	TrackNode *m_current;
};
