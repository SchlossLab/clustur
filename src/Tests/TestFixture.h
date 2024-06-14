//
// Created by gregj on 4/30/2024.
//

#ifndef TESTFIXTURE_H
#define TESTFIXTURE_H

class TestFixture {
public:
    virtual ~TestFixture() = default;
    virtual void Setup() = 0;
    virtual void TearDown() = 0;
};
#endif //TESTFIXTURE_H
