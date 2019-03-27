#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <string.h>

using namespace std;

ofstream out("out.dat", ios::out);

void writeTrack(const vector<int> &track, vector<int> outTrack, int last_i)
{
	for (unsigned i = last_i; i < track.size(); ++i)
	{
		int val = track.at(i);

		if (val == 0)
		{
			for (int num = 0; num < 6; ++num)
			{
				outTrack.push_back(num);
				writeTrack(track, outTrack, i+1);
				outTrack.pop_back();
			}

			return;
		}
		else if (val == 7)
		{
			for (int num = 12; num < 18; ++num)
			{
				outTrack.push_back(num);
				writeTrack(track, outTrack, i+1);
				outTrack.pop_back();
			}

			return;
		}
		else
		{
			outTrack.push_back(track.at(i)+5);
		}
	}

	out << outTrack.at(0);

	for (unsigned i = 1; i < outTrack.size(); ++i)
	{
		out << ' ' << outTrack.at(i);
	}

	out << endl;
}

int main()
{
	ifstream in("in.dat", ios::in);
	bool ok = in.is_open();

	const int bufSize = 1024;
	char *buff = (char*)malloc(sizeof(char) * bufSize);
	char *pch, *trash;

	while (!in.eof())
	{
		vector<int> track;
		in.getline(buff, bufSize);
		pch = strtok(buff, " ");

		while (pch != NULL)
		{
			int val = strtol(pch, &trash, 10);
			track.push_back(val);
			pch = strtok(NULL, " ");
		}

		vector<int> outTrack;
		writeTrack(track, outTrack, 0);
	}

	in.close();
	out.close();

	cout << "Completed";
	return 0;
}
