//      tagging_system.h
//      
//      Copyright 2011 Piotr Bentkowski <bentkowski.piotr@gmail.com>
//      
//      This program is free software; you can redistribute it and/or modify
//      it under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 2 of the License, or
//      (at your option) any later version.
//      
//      This program is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.
//      
//      You should have received a copy of the GNU General Public License
//      along with this program; if not, write to the Free Software
//      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//      MA 02110-1301, USA.



#ifndef _TAGGING_SYSTEM_
#define _TAGGING_SYSTEM_

#include <stdio.h>

/**
 * @class Tagging_system
 * 
 * @brief A singleton class which generates an unique tag for a cell.
 * 
 */
class Tagging_system {
private:
    Tagging_system();
    virtual ~Tagging_system();
    unsigned long int theTag;

public:
    static Tagging_system* getInstance();
    void setValue(unsigned long int v);
    void release();

    unsigned long int getTag();

private:
    static Tagging_system* s_pInstance;
};

#endif //_TAGGING_SYSTEM_
