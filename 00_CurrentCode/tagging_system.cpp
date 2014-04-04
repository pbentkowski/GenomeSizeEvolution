//      tagging_system.cpp
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


#include "tagging_system.h"

/**
 * @brief Data collecting method. Initialising singleton of tagging system.
 */
Tagging_system* Tagging_system::s_pInstance = NULL;

/**
 * @brief Data collecting method. Constructor
 */
Tagging_system::Tagging_system() {
}

/**
 * @brief Data collecting method. Destructor
 */
Tagging_system::~Tagging_system() {
}

/**
 * @brief Data collecting method. Gets instance of tagging system.
 */
Tagging_system* Tagging_system::getInstance() {
    if (NULL == s_pInstance) {
        s_pInstance = new Tagging_system();
    }
    return s_pInstance;
}

/**
 * @brief Data collecting method. Releases instance of tagging system.
 */
void Tagging_system::release() {
    if (NULL != s_pInstance) {
        delete s_pInstance;
        s_pInstance = NULL;
    }
}

/**
 * @brief Data collecting method. Sets a new value for tagging system.
 *
 * @param - number from which tagging begins.
 */
void Tagging_system::setValue(unsigned long v) {
    theTag = v;
}

/**
 * @brief Data collecting method. Obtaining an unique tag.
 *
 * @return - new unique tag (unsigned long int)
 */
unsigned long int Tagging_system::getTag() {
    theTag += 1;
    return theTag;
}
