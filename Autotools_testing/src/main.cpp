#include "hello.hpp"
#include "world.hpp"
#include <config.h> //Make configure results available
#include <stdio.h>
#include <iostream.h>

int main()
{

hello first_word;
world second_word;

}

std::cout<<PACKAGE_STRING; /* Use the preprocessor definitions from config.h */

first_word.print();
second_word.print();

return 0;

}

