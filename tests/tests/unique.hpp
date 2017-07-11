

#ifndef _UNIQUE_HPP_
#define _UNIQUE_HPP_

#include <iostream>
#include <vector>
#include <unordered_set>
#include <set>


#include "testslib.hpp"
#include "binarySearch.hpp"




    void unique(std::vector<double>& hashed_coords, std::unordered_set<double>& unique_hashes,
                std::vector<int>& unique_idx, std::vector<int>& idx)
    {
        unique_idx.clear();
        idx.clear();
        std::cout << "for 1" << std::endl;
        std::cout << "hashed_coords size" <<hashed_coords.size()<< std::endl;
        for (int i = 0; i < hashed_coords.size(); i++) {
            unique_hashes.insert(hashed_coords[i]);
        }
        unique_idx.resize(unique_hashes.size(),-1);
        // idx.resize(npixels);
        std::cout << "unique_hashes size" <<unique_hashes.size()<< std::endl;

        std::cout << "for 2" << std::endl;
        for (int i = 0; i < hashed_coords.size(); i++) {
            // int id = binarySearchRecursive(&unique_hashes[0],0,input.size(),hashed_coords[i]);
            std::unordered_set<double>::iterator got = unique_hashes.find (hashed_coords[i]);
            if(got != unique_hashes.end())
            {
                int id = std::distance(unique_hashes.begin(), got);
                idx.push_back(id);
                if(unique_idx[id] < 0) unique_idx[id] = i;
            }
        }

        std::cout << "for 2 end" << std::endl;

    }

    void unique(std::vector<double>& hashed_coords, std::vector<double>& unique_hashes,
                std::vector<int>& unique_idx,std::vector<int>& idx)
    {

        unique_idx.clear();
        idx.clear();
        unique_hashes.clear();


        std::set<double> input;
        std::cout << "for 1" << std::endl;
        std::cout << "hashed_coords size" <<hashed_coords.size()<< std::endl;
        for (int i = 0; i < hashed_coords.size(); i++) {
            // std::cout << "hashed_coords:"<<hashed_coords[i] << std::endl;
            input.insert(hashed_coords[i]);
        }
        unique_hashes.resize(input.size());
        unique_idx.resize(input.size(),-1);
        std::copy(input.begin(),input.end(),unique_hashes.begin());
        // std::cout << "input :" <<unique_hashes<< std::endl;
        std::cout << "input size" <<input.size()<< std::endl;

        std::cout << "for 2" << std::endl;
        for (int i = 0; i < hashed_coords.size(); i++) {
            int id = binarySearchRecursive(&unique_hashes[0],0,input.size(),hashed_coords[i]);
            if(id >= 0)
            {
                idx.push_back(id);
                if(unique_idx[id] < 0) unique_idx[id] = i;
            }
        }

        std::cout << "for 2 end" << std::endl;

    }

    void test_unique() {

        std::cout << "/ntest_unique/n" << std::endl;


        std::vector<double> hashed_coords = generateRandomVector<double>(npixels);
        // std::vector<double> hashed_coords(npixels,1.0);
        std::vector<double> unique_hashes;
        // std::unordered_set<double> unique_hashes;
        std::vector<int> unique_idx;
        std::vector<int> idx;
        unique(hashed_coords, unique_hashes, unique_idx, idx);

        std::cout << "hashed_coords:" << std::endl;
        PrintVector(hashed_coords);
        std::cout << "unique_hashes" << std::endl;
        PrintVector(unique_hashes);
        // PrintUnordered_set(unique_hashes);
        std::cout << "unique_idx" << std::endl;
        PrintVector(unique_idx);
        std::cout << "idx" << std::endl;
        PrintVector(idx);



    }






#endif //_UNIQUE_HPP_