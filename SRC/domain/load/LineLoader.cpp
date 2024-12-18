/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// Written: Brian Giffin & Himamshu Poudel, Oklahoma State University
//          09.2024
//
// Description: This file contains the class definition for LineLoader.

#include <LineLoader.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

Vector LineLoader::data(1);

LineLoader::LineLoader(int tag, int theElementTag)
  : ElementalLoad(tag, LOAD_TAG_LineLoader, theElementTag)
{
}

LineLoader::LineLoader()
  : ElementalLoad(LOAD_TAG_LineLoader)
{
}

LineLoader::~LineLoader()
{
}

const Vector &
LineLoader::getData(int &type, double loadFactor)
{
	type = LOAD_TAG_LineLoader;

	return data;
}

int
LineLoader::sendSelf(int commitTag, Channel &theChannel)
{
	int res = 0;
    static ID iddata(3);
    int dataTag = this->getDbTag();
    iddata(0) = this->getTag();
    iddata(1) = dataTag;
    iddata(2) = eleTag;

    res = theChannel.sendID(dataTag, commitTag, iddata);
    if (res < 0) {
        opserr << "WARNING LineLoader::sendSelf() - " << this->getTag() << " failed to send iddata\n";
        return res;
    }

    return res;
}

int
LineLoader::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int res = 0;
    static ID iddata(3);
    int dataTag = this->getDbTag();

    res = theChannel.recvID(dataTag, commitTag, iddata);
    if (res < 0) {
        opserr << "WARNING LineLoader::recvSelf() - " << this->getTag() << " failed to receive iddata\n";
        return res;
    }
     this->setTag(iddata(0));
     eleTag = iddata(2);

	return res;
}

void
LineLoader::Print(OPS_Stream &s, int flag)
{
	s << "LineLoader...";
	s << "  element acted on: " << eleTag << endln;
}

