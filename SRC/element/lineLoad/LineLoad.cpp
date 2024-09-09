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
                                                                        
                                                                        
// Written: Brian Giffin & Himamshu Poudel
// Created: 09.2024

// Description: This file contains the implementation for the LineLoad class.

#include "LineLoad.h"
#include <Information.h>
#include <ElementResponse.h>
#include <ID.h> 
#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <G3Globals.h>
#include <ErrorHandler.h>
#include <NDMaterial.h>
#include <ElementalLoad.h>
#include <packages.h>

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

Matrix LineLoad::tangentStiffness(LL_NUM_DOF, LL_NUM_DOF);
Vector LineLoad::internalForces(LL_NUM_DOF);
Vector LineLoad::theVector(LL_NUM_DOF);

#include <elementAPI.h>
static int num_LineLoad = 0;
static std::string globalLoadLibName = "";
static void* globalLoadLibPtr = nullptr;
static LineLoadFunct globalLoadFunctPtr = nullptr;

void *
OPS_LineLoad(void)
{
  if (num_LineLoad == 0) {
    num_LineLoad++;
    opserr<<"LineLoad element - Written: B.Giffin, H.Poudel, Oklahoma State University\n";
  }

  // Pointer to a uniaxial material that will be returned
  Element *theElement = 0;

  if (OPS_GetNumRemainingInputArgs() != 5) {
    opserr << "Want: element LineLoad eleTag? iNode? jNode? radius? lib?\n";
    return 0;
  }
    
  int    iData[3];
  double dData[1];
  const char* sData;

  int numData = 3;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element LineLoadElement" << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data: element LineLoad " << iData[0] << endln;
    return 0;	
  }

  sData = OPS_GetString();

  // Parsing was successful, allocate the material
  theElement = new LineLoad(iData[0], iData[1], iData[2], dData[0], sData);

  if (theElement == 0) {
    opserr << "WARNING could not create element of type LineLoadElement\n";
    return 0;
  }

  return theElement;
}

// constructors:
LineLoad::LineLoad(int tag, int Nd1, int Nd2, double radius, const char* lib)
 :Element(tag,ELE_TAG_LineLoad),     
   myExternalNodes(LL_NUM_NODE),
   dcrd1(LL_NUM_NDF),
   dcrd2(LL_NUM_NDF),
   libName(lib)
{
    myExternalNodes(0) = Nd1;
    myExternalNodes(1) = Nd2;

    my_radius = radius;

	mLoadFactor = 1.0;

    int ret = dynamicLibraryLoad();
    if (ret != 0) {
      opserr << "WARNING could not load dynamic library function for LineLoadElement\n";
    }
}

LineLoad::LineLoad()
  :Element(0,ELE_TAG_LineLoad),     
   	myExternalNodes(LL_NUM_NODE),
   	dcrd1(LL_NUM_NDF),
        dcrd2(LL_NUM_NDF),
        libName()
{
}

//  destructor:
LineLoad::~LineLoad()
{
}

int
LineLoad::getNumExternalNodes(void) const
{
    return LL_NUM_NODE;
}

const ID &
LineLoad::getExternalNodes(void) 
{
    return myExternalNodes;
}

Node **
LineLoad::getNodePtrs(void)
{
    return theNodes;                        
}

int
LineLoad::getNumDOF(void) 
{
    return LL_NUM_DOF;
}

void
LineLoad::setDomain(Domain *theDomain)
{
    theNodes[0] = theDomain->getNode(myExternalNodes(0));
    theNodes[1] = theDomain->getNode(myExternalNodes(1));

    for (int i = 0; i < LL_NUM_NODE; i++) {
    	if (theNodes[i] == 0)
        	return;  // don't go any further - otherwise segmentation fault
    }

    dcrd1 = theNodes[0]->getCrds();
    dcrd2 = theNodes[1]->getCrds();

    // call the base class method
    this->DomainComponent::setDomain(theDomain);
}

int
LineLoad::commitState()
{
	int retVal = 0;
    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
    	opserr << "LineLoad::commitState () - failed in base class";
    }    

    return 0; 
}

int
LineLoad::revertToLastCommit()
{
	return 0;
}

int
LineLoad::revertToStart()
{
	return 0;
}

int
LineLoad::update(void)
{
	return 0;
}

int
LineLoad::dynamicLibraryLoad()
// this function loads the dynamic library function for applying line loads to the member
{
 
    // try existing loaded routine
    if (libName == globalLoadLibName) {

      // set the library handle and function ptr &  return it
      libHandle = globalLoadLibPtr;
      lineLoadFunctPtr = globalLoadFunctPtr;
      return 0;
      
    }

    // ty to load new routine from dynamic library in load path
    const char *libraryName  = libName.c_str();
    const char *functionName = "OPS_ApplyLineLoad";

    int res = getLibraryFunction(libraryName, functionName, &libHandle, (void**)&lineLoadFunctPtr);

    if (res == 0) {
      // set the global data consistent with the newly loaded library routine
      globalLoadLibName = libName;
      globalLoadLibPtr = libHandle;
      globalLoadFunctPtr = lineLoadFunctPtr;
      return 0;
    } else {
      return 1;
    }
    
}

const Matrix &
LineLoad::getTangentStiff(void)
{
    tangentStiffness.Zero();
    return tangentStiffness;
}

const Matrix &
LineLoad::getInitialStiff(void)
{
    return getTangentStiff();
}
    
void 
LineLoad::zeroLoad(void)
{
    return;
}

int 
LineLoad::addLoad(ElementalLoad *theLoad, double loadFactor)
{
	int type;
	const Vector &data = theLoad->getData(type, loadFactor);

	if (type == LOAD_TAG_LineLoader) {
		mLoadFactor = loadFactor;
		return 0;
	} else {
		opserr << "LineLoad::addLoad() - ele with tag: " << this->getTag() << " does not accept load type: " << type << endln;
		return -1;
	}

	return -1;
}

int 
LineLoad::addInertiaLoadToUnbalance(const Vector &accel)
{
	return 0;
}

const Vector &
LineLoad::getResistingForce()
{
	internalForces.Zero();

	// get current time
	Domain *theDomain = this->getDomain();
	double t = theDomain->getCurrentTime();

	// initialize the arrays of nodal forces and coordinates
	double forces[LL_NUM_NODE*LL_NUM_NDF] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	double coordinates[LL_NUM_NODE*LL_NUM_NDF] = { dcrd1(0), dcrd1(1), dcrd1(2), dcrd2(0), dcrd2(1), dcrd2(2) };
	
	// call the external library function to determine the resulting nodal forces applied to the element
	lineLoadFunctPtr(&forces[0], &coordinates[0], this->getTag(), my_radius, t);
	
	// compute internal forces by scaling the nodal forces by the associated load factor
	
	// loop over nodes
	for(int j = 0; j < LL_NUM_NODE; j++) {
	  // loop over dof
	  for(int k = 0; k < LL_NUM_NDF; k++) {
	    internalForces[j*LL_NUM_NDF+k] = internalForces[j*LL_NUM_NDF+k] - mLoadFactor*forces[j*LL_NUM_NDF+k];
	  }
	}

	return internalForces;
}

const Vector &
LineLoad::getResistingForceIncInertia()
{       
  	return getResistingForce();
}

int
LineLoad::sendSelf(int commitTag, Channel &theChannel)
{
  int res;

  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();

  // LineLoad packs its data into a Vector and sends this to theChannel
  // along with its dbTag and the commitTag passed in the arguments

  static Vector data(3 + 2*LL_NUM_NDF);
  data(0) = this->getTag();
  data(1) = my_radius;
  data(2) = mLoadFactor;

  for (int i = 0; i < LL_NUM_NDF; i++) {
    data(3+0*LL_NUM_NDF+i) = dcrd1(i);
    data(3+1*LL_NUM_NDF+i) = dcrd2(i);       
  }
  
  res = theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr <<"WARNING LineLoad::sendSelf() - " << this->getTag() << " failed to send data\n";
    return -1;
  }           

  // LineLoad then sends the tags of its four nodes
  res = theChannel.sendID(dataTag, commitTag, myExternalNodes);
  if (res < 0) {
    opserr <<"WARNING LineLoad::sendSelf() - " << this->getTag() << " failed to send myExternalNodes\n";
    return -2;
  }

  return 0;
}

int
LineLoad::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res;
  int dataTag = this->getDbTag();

  // LineLoad creates a Vector, receives the Vector and then sets the 
  // internal data with the data in the Vector
  static Vector data(3 + 2*LL_NUM_NDF );
  res = theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr <<"WARNING LineLoad::recvSelf() - failed to receive Vector\n";
    return -1;
  }           
  
  this->setTag(int(data(0)));
  my_radius = data(1);
  mLoadFactor = data(2);

  for (int i = 0; i < LL_NUM_NDF; i++) {
    dcrd1(i)  = data(3+0*LL_NUM_NDF+i);
    dcrd2(i)  = data(3+1*LL_NUM_NDF+i);
  }

  // LineLoad now receives the tags of its four external nodes
  res = theChannel.recvID(dataTag, commitTag, myExternalNodes);
  if (res < 0) {
    opserr <<"WARNING LineLoad::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return -2;
  }

  return 0;
}

int
LineLoad::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
  return 0;
}

void
LineLoad::Print(OPS_Stream &s, int flag)
{
    opserr << "LineLoad, element id:  " << this->getTag() << "\n";
    opserr << "   Connected external nodes:  " ; 
    for (int i = 0; i<LL_NUM_NODE; i++) {
    	opserr << myExternalNodes(i)<< " ";
    }
    opserr << "\n";
    return;
}

Response*
LineLoad::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  return Element::setResponse(argv, argc, output);
}

int 
LineLoad::getResponse(int responseID, Information &eleInfo)
{
  return Element::getResponse(responseID, eleInfo);
}

