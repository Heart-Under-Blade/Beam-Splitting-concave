#include "Tracks.h"

#include <iostream>

#include "Beam.h"

int Tracks::FindGroupByTrackId(const IdType &trackId) const
{
	for (size_t i = 0; i < size(); ++i)
	{
		for (int j = 0; j < (*this)[i].size; ++j)
		{
			if ((*this)[i].arr[j] == trackId)
			{
				return (*this)[i].groupID;
			}
		}
	}

	if (size() == 0)
	{
		return 0;
	}

	return -1;
}

void Tracks::ImportTracks(int nFacets, const std::string &filename)
{
	const int bufSize = 1024;
	std::ifstream trackFile(filename, std::ios::in);

	if (!trackFile.is_open())
	{
		std::cerr << "Track file not found" << std::endl;
		throw std::exception();
	}

	char *buff = (char*)malloc(sizeof(char) * bufSize);

	TrackGroup buffGroup;

	while (!trackFile.eof())
	{
		trackFile.getline(buff, bufSize);

		vector<int> track;

		char *ptr, *trash;
		ptr = strtok(buff, " ");

		size_t groupIndex = 0;
		bool haveGroup = false;

		while (ptr != NULL)
		{
			if (ptr[0] == ':')
			{
				haveGroup = true;
				ptr = strtok(NULL, " ");
				groupIndex = strtol(ptr, &ptr, 10);

				if (groupIndex >= size())
				{
					for (size_t i = size(); i <= groupIndex; ++i)
					{
						(*this).push_back(TrackGroup());
					}
				}

				(*this)[groupIndex].groupID = groupIndex;
				break;
			}

			int tmp = strtol(ptr, &trash, 10);
			track.push_back(tmp);
			ptr = strtok(NULL, " ");
		}

		int trackID = 0;

		for (int t : track)
		{
			trackID += (t + 1);
			trackID *= (nFacets + 1);
		}

		if (haveGroup)
		{
			(*this)[groupIndex].tracks.push_back(track);
			(*this)[groupIndex].arr[(*this)[groupIndex].size++] = trackID;
		}
		else
		{
			buffGroup.tracks.push_back(track);
			buffGroup.arr[buffGroup.size++] = trackID;
		}

		track.clear();
	}

	if (buffGroup.size != 0) // добавляем треки без группы в отдельные группы
	{
		for (int i = 0; i < buffGroup.size; ++i)
		{
			TrackGroup newGroup;
			newGroup.arr[newGroup.size++] = buffGroup.arr[i];
			newGroup.tracks.push_back(buffGroup.tracks[i]);
			newGroup.groupID = size();
			push_back(newGroup);
		}
	}
}

void Tracks::RecoverTrack(const Beam &beam, int facetNum,
						  std::vector<int> &track)
{
	int coef = facetNum + 1;
	std::vector<int> tmp_track;

	auto tmpId = beam.id/coef;

	for (int i = 0; i <= beam.nActs; ++i)
	{
#ifdef _DEBUG // DEB
		int tmp = (tmpId%coef);
#else
		int tmp = (tmpId%coef).toInt();
#endif
		tmpId -= tmp;
		tmpId /= coef;
		tmp -= 1;
		tmp_track.push_back(tmp);
	}

	for (int i = tmp_track.size()-1; i >= 0; --i)
	{
		track.push_back(tmp_track.at(i));
	}
}
