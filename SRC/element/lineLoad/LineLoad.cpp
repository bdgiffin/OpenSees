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

#include <math.h>
#include <stdlib.h>
#include <stdio.h> 

double LineLoad :: oneOverRoot3 = 1.0/sqrt(3.0);
double LineLoad :: GsPts[2];

Matrix LineLoad::tangentStiffness(LL_NUM_DOF, LL_NUM_DOF);
Vector LineLoad::internalForces(LL_NUM_DOF);
Vector LineLoad::theVector(LL_NUM_DOF);

#include <elementAPI.h>
static int num_LineLoad = 0;

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
   g1(LL_NUM_NDF),
   myNI(LL_NUM_NODE),
   dcrd1(LL_NUM_NDF),
   dcrd2(LL_NUM_NDF)
{
    myExternalNodes(0) = Nd1;
    myExternalNodes(1) = Nd2;

	GsPts[0] = -oneOverRoot3;
	GsPts[1] = +oneOverRoot3;

    my_radius = radius;

	mLoadFactor = 1.0;
}

LineLoad::LineLoad()
  :Element(0,ELE_TAG_LineLoad),     
   	myExternalNodes(LL_NUM_NODE),
   	g1(LL_NUM_NDF),
   	myNI(LL_NUM_NODE),
   	dcrd1(LL_NUM_NDF),
   	dcrd2(LL_NUM_NDF)
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
LineLoad::UpdateBase(double Xi)
// this function calculates g1 and NI for given Xi
{
    double oneMinusXi  = 1 - Xi;
    double onePlusXi   = 1 + Xi;

    // calculate vector g1
    // g1 = d(x_Xi)/dXi
    g1 = (dcrd2 - dcrd1) * 0.5;

	// shape functions
	myNI(0) = 0.5 * oneMinusXi;
	myNI(1) = 0.5 * onePlusXi;

    return 0;
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

	// loop over Gauss points
	for(int i = 0; i < 2; i++) {
		this->UpdateBase(GsPts[i]);

		// get the current position along the length of the element
		Vector icrd = myNI(0)*dcrd1 + myNI(1)*dcrd2;

		// -------------------- EXAMPLE DRAG FORCE CALCULATION -------------------- //

		// set the ambient velocity and density at the current location
		double drag_coeff = 1.0;
		double area = 2.0*my_radius;
		double density = 0.001;
		Vector velocity(3);
		double v_ref = 1.0;
		double time_scaling = 1.0 - exp(-t);
		velocity(0) = v_ref*time_scaling*(-icrd(1));
		velocity(1) = v_ref*time_scaling*(+icrd(0));
		velocity(2) = v_ref*time_scaling*(1.0 - exp(-icrd(2)));

		// get the current element length and orientation vector
		Vector lambda = 2.0*g1;
		double length = sqrt(lambda(0)*lambda(0) + lambda(1)*lambda(1) + lambda(2)*lambda(2));
		lambda = lambda/length;

		// get the projected wind velocity, removing the axial component
		Vector proj_velocity = velocity - lambda*(lambda(0)*velocity(0) + lambda(1)*velocity(1) + lambda(2)*velocity(2));
		double norm_proj_vel = sqrt(proj_velocity(0)*proj_velocity(0) + proj_velocity(1)*proj_velocity(1) + proj_velocity(2)*proj_velocity(2));

		// compute the total drag load
		Vector drag_force = (0.5*density*drag_coeff*area*length*area*norm_proj_vel)*proj_velocity;
		
		// ------------------------------------------------------------------------ //

		// loop over nodes
		for(int j = 0; j < LL_NUM_NODE; j++) {
			// loop over dof
			for(int k = 0; k < LL_NUM_NDF; k++) {
			        internalForces[j*LL_NUM_NDF+k] = internalForces[j*LL_NUM_NDF+k] - mLoadFactor*myNI(j)*drag_force(k);
			}
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

  static Vector data(3 + 3*LL_NUM_NDF + LL_NUM_NODE);
  data(0) = this->getTag();
  data(1) = my_radius;
  data(2) = mLoadFactor;

  for (int i = 0; i < LL_NUM_NDF; i++) {
    data(3+             i) = g1(i);
    data(3+1*LL_NUM_NDF+i) = dcrd1(i);
    data(3+2*LL_NUM_NDF+i) = dcrd2(i);       
  }
  for (int i = 0; i < LL_NUM_NODE; i++)
    data(3+3*LL_NUM_NDF+i) = myNI(i);
  
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
  static Vector data(3 + 3*LL_NUM_NDF + LL_NUM_NODE);
  res = theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr <<"WARNING LineLoad::recvSelf() - failed to receive Vector\n";
    return -1;
  }           
  
  this->setTag(int(data(0)));
  my_radius = data(1);
  mLoadFactor = data(2);

  for (int i = 0; i < LL_NUM_NDF; i++) {
    g1(i)     = data(3+             i);
    dcrd1(i)  = data(3+1*LL_NUM_NDF+i);
    dcrd2(i)  = data(3+2*LL_NUM_NDF+i);
  }
  for (int i = 0; i < LL_NUM_NODE; i++)
    myNI(i) = data(3+3*LL_NUM_NDF+i);

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

