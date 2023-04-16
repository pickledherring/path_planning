#include <iostream>
#include <vector>
using namespace std;
 
int* fun()
{
    int* arr = new int[2];
 
    /* Some operations on arr[] */
    arr[0] = 10;
    arr[1] = 20;
 
    return arr;
}
 
int main()
{
    // int* ptr = fun();
    // cout << ptr[0] << ", " << ptr[1]<<endl;
    // delete[] ptr;
    vector<int*> vec;
    for (int i = 0; i < 3; i++) {
        int* pos = new int[2];
        pos[0] = i;
        pos[1] = i;
        vec.push_back(pos);
    }
    int new_pos[2];
    vec.pop_back();
    for (int i = 0; i < vec.size(); i ++) {
        cout<<"vec ["<<i<<"]: "<<vec[i][0]<<", "<<vec[i][1]<<endl;
    }
    vec.pop_back();
    new_pos[0] = vec.back()[0];
    new_pos[1] = vec.back()[1];
    
    cout<<new_pos[0]<<", "<<new_pos[1];

    return 0;
}