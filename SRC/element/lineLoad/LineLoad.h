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
                                                                        
#ifndef LineLoad_h
#define LineLoad_h

// Written: Brian Giffin & Himamshu Poudel
// Created: 09.2024

// Description: This file contains the class definition for LineLoad. 

#include <Element.h>
#include <Node.h>
#include <Vector.h>
#include <Matrix.h>
#include <NDMaterial.h>
#include <ID.h>

#include <string>

class Domain;
class Node;
class Channel;
class NDMaterial;
class FEM_ObjectBroker;

typedef void (*LineLoadFunct)(double* forces, const double* coordinates, int element_tag, double radius, double time);

class LineLoad : public Element
{
  public:
    LineLoad(int tag, int Nd1, int Nd2, double radius, const char* lib); 
    LineLoad();
    ~LineLoad();

    // public methods to obtain information about dof & connectivity    
    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);
    int getNumDOF(void);	
    void setDomain(Domain *theDomain);

    // public methods to set the state of the element    
    int commitState(void);
    int revertToLastCommit(void);        
    int revertToStart(void);        
    int update(void);
    
    // public methods to obtain stiffness, mass, damping and 
    // residual information    
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);    

    void zeroLoad(void);	
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);
    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);            

    // public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);
    int displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode);

    void Print(OPS_Stream &s, int flag =0);    

    Response *setResponse(const char **argv, int argc, OPS_Stream &output);
    int getResponse(int responseID, Information &eleInformation);

  protected:
    
  private:

  enum {LL_NUM_NODE = 2}; // number of nodes per element
  enum {LL_NUM_NDF = 3}; // d.o.f. per node
  enum {LL_NUM_DOF = 6}; // degrees of freedom per element
  enum {LL_NUM_DDOF = 6}; // displacement degrees of freedom per element
  
    // method to load the dynamic library line load routine
    int dynamicLibraryLoad();

    ID  myExternalNodes;      // contains the tags of the end nodes
    static Matrix tangentStiffness;  // Tangent Stiffness matrix
    static Vector internalForces;    // vector of Internal Forces
    static Vector theVector;         // vector to return the residual

    double my_radius;       // the effective radius of the structural element

    Node *theNodes[LL_NUM_NODE];

    Vector dcrd1;             // current coordinates of node 1
    Vector dcrd2;             // current coordinates of node 2

    std::string libName;       // the user-defined external library name
    void* libHandle;
    LineLoadFunct lineLoadFunctPtr;

	double mLoadFactor;       // factor from load pattern
};

#endif




