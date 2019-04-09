#pragma once

#include "Facet.h"

#include <vector>

struct TrackNode
{
	bool isCompleted;
	Facet *m_facet;
	TrackNode* m_parent;
	std::vector<TrackNode*> children;

	TrackNode(Facet *facet, TrackNode* parent = nullptr)
	{
		isCompleted = false;
		m_facet = facet;
		m_parent = parent;
	}

	~TrackNode()
	{
		for (auto child : children)
		{
			delete child;
		}
	}

	TrackNode *AddChild(Facet *facet)
	{
		TrackNode* child = new TrackNode(facet, this);
		children.push_back(child);
		return child;
	}

	TrackNode *FindNode(int id)
	{
		for (auto child : children)
		{
			if (child->m_facet->index == id)
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
