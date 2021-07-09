
#include <iostream>
#include <vector>
#include <cctype>

using namespace std;

pair<vector<char>, vector<char>> ReadTexts() {

	int textSize;
	char a;
	pair<vector<char>, vector<char>> t;
	cin >> textSize;

	
	for(int i = 0; i < textSize; ++i) {
		a = getchar();
		t.first.push_back(tolower(a));
    }
	a = getchar();
	for(int i = 0; i < textSize; ++i) {
		a = getchar();
		t.second.push_back(tolower(a));
	}

	
	return t;
}

double TextCompare(pair<vector<char>, vector<char>> &texts) {
	int count = 0;
	for (int i = 0; i < texts.first.size(); i++) {
		if (texts.first[i] == texts.second[i]) {
			count++;
		}
	}
	return (double)count / texts.first.size();
}

int main() {
	
	pair<vector<char>, vector<char>> texts = ReadTexts();
	cout << TextCompare(texts);

	return 0;
}