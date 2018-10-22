/*************************************************************************
POMATO: POincare MAp TOpology: Extraction and visualization of the
        topological structure of the circular restricted three-body 
        problem for orbital mechanics applications. 

Authors: Wayne Schlei and Xavier Tricoche 

Copyright (c) 2013-2018, Purdue University 

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
**************************************************************************/


//Avizo
#include <QApplication>
// #include <hxcore/HxMessage.h>
// #include <hxcore/HxProgressInterface.h>
// #include <hxcore/HxController.h>
// #include <hxcore/HxObjectPool.h>
// #include <hxfield/HxCoordType.h>
#include <mypackage/QStringInterface.h>
#include <mypackage/CR3BPsystem.h>
// #include <mypackage/HxTrajectory.h>
// #include <mypackage/HxInitialState.h>
// #include <mypackage/HxManeuver.h>
#include <poincare/MyLine.h>
#include <poincare/FixedPointData.h>
#include <poincare/ManifoldData.h>
#include <poincare/DisplayManifolds.h>

//Rendering Elements
// #include <Inventor/nodes/SoTransform.h>
// #include <Inventor/nodes/SoTranslation.h>
// #include <Inventor/SbRotation.h>
// #include <Inventor/nodes/SoSeparator.h>
// #include <Inventor/nodes/SoSphere.h>
// #include <Inventor/nodes/SoPointSet.h>
// #include <Inventor/nodes/SoLineSet.h>
// #include <Inventor/nodes/SoCoordinate3.h>
// #include <Inventor/nodes/SoMaterial.h>
// #include <Inventor/nodes/SoMaterialBinding.h>
// #include <Inventor/nodes/SoVertexProperty.h>
// #include <Inventor/nodes/SoDrawStyle.h>
// #include <Inventor/nodes/SoFont.h>
// #include <Inventor/nodes/SoText2.h>
// #include <Inventor/nodes/SoDepthOffset.h>
// #include <Inventor/nodes/SoMarkerSet.h>
// #include <Inventor/nodes/SoShapeHints.h>

//Events
// #include <Inventor/SoPickedPoint.h>
// #include <Inventor/events/SoLocation2Event.h>
// #include <Inventor/events/SoMouseButtonEvent.h>
// #include <Inventor/events/SoMouseWheelEvent.h>
// #include <Inventor/events/SoKeyboardEvent.h>
// #include <Inventor/details/SoLineDetail.h>
// #include <Inventor/actions/SoCallbackAction.h>


//Standard Libs
#include <vector>
#include <list>
#include <map>
#include <queue>
#include <boost/format.hpp>
#include <boost/limits.hpp>

//API - chris
#include <math/rational.hpp>
#include <math/fixed_vector.hpp>
#include <util/wall_timer.hpp>
//API - xavier
#include <maps/definitions.hpp>
#include <maps/map_analysis.hpp>
#include <maps/fixpoints.hpp>
#include <maps/poincare_map.hpp>
//API - orbital
#include <maps/DP45wrapper.hpp>
#include <cr3bp/cr3bp.hpp>
#include <cr3bp/tracker.hpp>
#include <cr3bp/planar_section.hpp>
#include <orbital/corrections.hpp>
#include <topology/invariant_manifold.hpp>
//#include <topology/ManifoldDataStorage.hpp>

#if _OPENMP
#include <omp.h>
#endif
using namespace xavier;
using namespace topology;
//using namespace std; //Some functions have titles the same as other std::functions();

//Type Definitions
typedef xavier::dp5wrapper<double, 42>                                ode_solver;
typedef orbital::cr3bp                                                rhs_type;
typedef orbital::planar_section<rhs_type, 6, 42>                      section_type;
typedef xavier::poincare_map<rhs_type, ode_solver, section_type >     map_type;
typedef xavier::state_info<rhs_type::value_type, rhs_type::dimension> return_state;
typedef map_type::return_type                                         return_type;
typedef nvis::fixed_vector<double,6>             vec6;
typedef nvis::fixed_vector<double,42>            vec42;
typedef nvis::fixed_matrix<double,6>             mat6;
typedef std::vector<nvis::vec2>                  orbit_type;
typedef std::vector<fp_chain>                    chain_type;
//typedef std::vector<Separatrix>                  sep_type;
typedef ManifoldData::MapDiscont                 MapDiscont;
typedef ManifoldData::ManifoldSeg                ManifoldSeg;



///// Callback when motion is detected
//void DisplayManifolds_clickMotionCB(void *userData, SoEventCallback *eventCallback)
//{
//    //Call this module
//    DisplayManifolds* module = (DisplayManifolds*) userData;
//    //Interface for QStrings
//    QStringInterface qi;
//    
//    //States of module
//    bool hasInput = module->hasInput;
//    bool buttonPressed = module->buttonPressed;
//    bool wwDesignStarted = module->wwTopoDesignON;
//    
//    //Event and state of events
//    const SoEvent *event = eventCallback->getEvent();
//    bool isButtonEvent = (event->isOfType(SoMouseButtonEvent::getClassTypeId()));
//    bool isLoc2Event = (event->isOfType(SoLocation2Event::getClassTypeId()));
//    
//    //Motion Events:
//    if ( isLoc2Event ) { 
//      SoLocation2Event *loc2Event = (SoLocation2Event*) event;
//      if(loc2Event->getEventSource() == SoLocation2Event::MOUSE_MOVE) { //Motion only
//        //Picking is handled automatically so no need to check motion event
//        const SoPickedPoint * pp = eventCallback->getPickedPoint();
//        if (pp!=NULL) { //If mouse intersects something
//            //Need a check to see if this contains the intended node in the path:
//            const SoPath *pPath = pp->getPath();
//            
//            //Looking for only the main manifold render node (may not be necessary!)
//            bool containsTarget = pPath->containsNode( module->targetNode );
//            /* if (containsTarget) {
//              theMsg->printf("Target Node (SoLineSet making scene) FOUND!");
//            } else 
//              theMsg->printf("Target Node LOST!");
//            */
//            
//            //Work through the items in the path to find the first line 
//            SoLineSet *theIntersectedLine = 0;
//            for(int i=0; i<pPath->getLength(); i++) {
//              SoNode *currentNode = pPath->getNode(i);
//              //Stop at first lineset
//              if (currentNode->getTypeId() == SoLineSet::getClassTypeId()) {
//                theIntersectedLine = (SoLineSet*) currentNode;
//                break;
//              }
//            }
//            
//            // -------------------------------------------------------------
//            // Standard Interaction: Hover over manifold point WITHOUT click
//            //  + Heteroclinic option
//            //  + Topology design from Initial state
//            // -------------------------------------------------------------
//            //If it intersects lines of this scene's lineset
//            if ( !buttonPressed && containsTarget && (theIntersectedLine != NULL) ) {
//              //Extract information from picked point and SoLineSet
//              SoLineDetail *lineDeets = (SoLineDetail*) pp->getDetail(); //Detail about hit node
//              int32_t lineIdx = lineDeets->getLineIndex(); //tail 
//              int32_t partIdx = lineDeets->getPartIndex(); //tail
//              SbVec3f intersectionPoint = pp->getPoint();
//              //Need to test overall ray pick action path to see if it hits both manifold objects
//              const SoPickedPointList ppList = eventCallback->getAction()->getPickedPointList();
//              bool containsUnstable = false; int uPathIdx = 0;
//              bool containsStable = false; int sPathIdx = 0;
//              bool containsDeltaVLine = false;
//              //Loop through list to see if it holds both nodes
//              for(int i=0; i<ppList.getLength(); i++) {
//                if(!containsUnstable) {
//                  containsUnstable = ppList[i]->getPath()->containsNode( module->wuTarget );
//                  if (containsUnstable) {
//                    uPathIdx = i;
//                    //theMsg->printf("UnstableManifold Node found: index = %d", uPathIdx);
//                  }
//                }
//                if(!containsStable) {
//                  containsStable = ppList[i]->getPath()->containsNode( module->wsTarget );
//                  if (containsStable) {
//                    sPathIdx = i;
//                    //theMsg->printf("StableManifold Node found: index = %d", sPathIdx);
//                  }
//                }
//                //Also need to look for DeltaVLine node for topology design platform
//                if(hasInput && !containsDeltaVLine) {
//                  containsDeltaVLine = ppList[i]->getPath()->containsNode( module->topoDesignTargetNode);
//                }
//                
//              }
//              //Check for Heteroclinic/Homoclinic connections by seeing if both are detected
//              bool hcDetected = false;
//              if( !hasInput && containsUnstable && containsStable ) {
//                //Skip if we have InitialState input for Topology Design platform
//                //theMsg->printf("Possible Heteroclinic/homoclinic connection located");
//                //getPickedPoint() only grabs the first node it hits, so we have to get details another way
//                //We have to determine what is what of this selection!
//                SoLineDetail *ulineDeets = (SoLineDetail*) ppList[uPathIdx]->getDetail();
//                SoLineDetail *slineDeets = (SoLineDetail*) ppList[sPathIdx]->getDetail();
//                //Force selection to the unstable manifold 
//                module->selectedUnstable = true;
//                //Gather end point information within module (send line/part ids)
//                module->selectedLineID = (int) ulineDeets->getLineIndex();
//                module->selectedPartID = (int) ulineDeets->getPartIndex();
//                //Have to get the partner node with Line and Part 
//                module->hcPartnerLineID = (int) slineDeets->getLineIndex();
//                module->hcPartnerPartID = (int) slineDeets->getPartIndex();
//                //printf("At HC selection: slID = %d spID = %d   | hclID = %d hcpID = %d\n",
//                //        module->selectedLineID, module->selectedPartID, 
//                //        module->hcPartnerLineID, module->hcPartnerPartID);
//                //Only if we found a valid connection (segs intersect)
//                if (module->makeHHConnection()) { 
//                  //Indicate that we have detected a connection
//                  hcDetected = true;
//                  module->hcOutput.setValue(0,1); //This calls compute()
//                }
//                //Pass the picked point to module
//                module->pickedPoint = intersectionPoint;
//                module->show = true;
//              }
//              //Reset the picked point and path information (as things were likely shifted)
//              // ->This does not apply to Topology Design except for prior to a selection
//              const SoPath *pPath2 = eventCallback->getPickedPoint()->getPath();
//              if ( !hasInput && !hcDetected && pPath2->containsNode( module->wuTarget ) ) {
//                module->selectedUnstable = true;
//                module->hcOutput.setValue(0,0);
//                module->hcPartnerLineID = 0;
//                module->hcPartnerPartID = 0;
//                //printf(" hcValues cleared LINE %d\n",__LINE__);
//                //Gather end point information within module (send line/part ids)
//                module->selectedLineID = (int) lineIdx;
//                module->selectedPartID = (int) partIdx;
//                //Pass the picked point to module
//                module->pickedPoint = intersectionPoint;
//                module->show = true;
//              } else if ( !hasInput && !hcDetected && pPath2->containsNode( module->wsTarget ) ) {
//                module->selectedUnstable = false;
//                module->hcOutput.setValue(0,0);
//                module->hcPartnerLineID = 0;
//                module->hcPartnerPartID = 0;
//                //printf(" hcValues cleared LINE %d\n",__LINE__);
//                //Gather end point information within module (send line/part ids)
//                module->selectedLineID = (int) lineIdx;
//                module->selectedPartID = (int) partIdx;
//                //Pass the picked point to module
//                module->pickedPoint = intersectionPoint;
//                module->show = true;
//              }
//              //std::cout << " LineIndex = " << (int) lineIdx << "  Part Index = " << (int) partIdx << "\n";
//              // -> Need lineIndex : which line within SoLineSet
//              QString lidStr; lidStr.setNum( lineIdx );
//              // -> Need partIndex : which segment within SoLineSet
//              QString partStr; partStr.setNum( partIdx );
//              QString manTypeStr((module->selectedUnstable)? "(Wu)": "(Ws)"); 
//              QString hcStr((module->hcOutput.getValue(0)==1)? " (HC)": "" );
//              QString infoStr =  manTypeStr + " LineIndex = " + lidStr + " PartIndex = " + partStr + hcStr;
//              module->selectionInfo.setValue( qi.setString(infoStr) );
//              
//              //Topology Design options forces select choices
//              bool wuDesign = false, wsDesign = false;
//              //Topology design version of selection must access the correct manifold
//              if (hasInput) {
//                //Always pass the picked point to the module
//                module->pickedPoint = intersectionPoint;
//                switch(module->portDesignOpt.getValue()) {
//                  case 0:
//                    //Looking for stable manifold (forward time)
//                    wsDesign = true; wuDesign = false;
//                    break;
//                  case 3 :
//                    //Looking for stable manifold (backward time)
//                    wsDesign = true; wuDesign = false;
//                    break;
//                  default :
//                    //Assume others are looking for unstable
//                    wuDesign = true; wsDesign = false;
//                }
//                if(wuDesign && containsDeltaVLine && containsUnstable) {
//                  //UNSTABLE manifold needed that intersects design line:
//                  //We have to determine what is what of this selection!
//                  SoLineDetail *ulineDeets = (SoLineDetail*) ppList[uPathIdx]->getDetail();
//                  //Force selection to the unstable manifold 
//                  module->selectedUnstable = true;
//                  //Gather end point information within module (send line/part ids)
//                  module->selectedLineID = (int) ulineDeets->getLineIndex();
//                  module->selectedPartID = (int) ulineDeets->getPartIndex();
//                  //Compute the maneuver
//                  if( module->topologyDeltaV() ) {
//                      module->topoDesignSelected = true;
//                  } else {
//                      module->topoDesignSelected = false;
//                  }
//                }
//                else if (wsDesign && containsDeltaVLine && containsStable) {
//                  //STABLE manifold needed that intersects design line:
//                  SoLineDetail *slineDeets = (SoLineDetail*) ppList[sPathIdx]->getDetail();
//                  //Force selection to the stable manifold
//                  module->selectedUnstable = false;
//                  //Have to get the partner node with Line and Part 
//                  module->selectedLineID = (int) slineDeets->getLineIndex();
//                  module->selectedPartID = (int) slineDeets->getPartIndex();
//                  //Compute the maneuver
//                  if( module->topologyDeltaV() ) {
//                      module->topoDesignSelected = true;
//                  } else {
//                      module->topoDesignSelected = false;
//                  }
//                } else {
//                  //Indicate that we don't yet have a valid design point 
//                  module->topoDesignSelected = false;
//                }
//                //Update the topology Design overlay rendering (always based on motion)
//                module->renderTopoDesignScene();
//                //Disable Heteroclinic/homoclinic operations
//                module->hcPartnerLineID = 0;
//                module->hcPartnerPartID = 0;
//                module->hcOutput.setValue(0,0);
//                //printf(" hcValues cleared LINE %d\n",__LINE__);
//              } //End Topology Design from Input
//
//              //Enforce simulation to be called
//              module->integrate = true;
//            } 
//            //else module->show = false; //Not the right stuff to show!
//          
//            // --------------------------------------------------------------------------------
//            // Motion Events for WtoW design (WITH button press)
//            // --------------------------------------------------------------------------------
//            if (buttonPressed && !hasInput) {
//              //IF DESIGN IS STARTED ---------------------------------------------------
//              if(wwDesignStarted) { //And drag is applied through button press
//                  //Get the picked point
//                  SbVec3f thePoint = pp->getPoint();
//                  //Always send to scene for rendering
//                  module->pickedPoint = thePoint;
//                  
//                  //To continue design, Look for intersecting objects of relevance
//                  const SoPickedPointList ppList = eventCallback->getAction()->getPickedPointList();
//                  bool wantUnstable = false;
//                  if(module->portDesignOpt.getValue() == 5) wantUnstable = true;
//                  bool containsWTarget = false; int wPathIdx = 0;
//                  bool containsDeltaVLine = false; 
//                  for(int i=0;i<ppList.getLength();i++) {
//                      //Manifold
//                      if(!containsWTarget) {
//                        if (wantUnstable) {
//                          //Look for unstable manifold node
//                          containsWTarget = ppList[i]->getPath()->containsNode( module->wuTarget );
//                        } else {
//                          //Look for stable manifold node
//                          containsWTarget = ppList[i]->getPath()->containsNode( module->wsTarget );
//                        }
//                        if (containsWTarget) wPathIdx = i;
//                      }
//                      //DeltaV constraint line
//                      if(!containsDeltaVLine) {
//                        containsDeltaVLine = ppList[i]->getPath()->containsNode( module->topoDesignTargetNode );
//                      }
//                  }
//                  //Debug:
//                  //theMsg->printf("PickedPoint intersection:  DVLine = %d  TargetManifold (%s) = %d",
//                  //                        (containsDeltaVLine)?1:0, (wantUnstable)?"Wu":"Ws", (containsWTarget)?1:0);
//                  
//                  //Must intersect both topoDesignTargetNode and the opposing manifold
//                  // to create and store a design
//                  if (containsDeltaVLine && containsWTarget) {
//                      //Tell module to compute design with given point node (temporarily)
//                      SoLineDetail *lineDeets = (SoLineDetail*) ppList[wPathIdx]->getDetail();
//                      module->selectedUnstable = wantUnstable;
//                      //Gather point information from line details (for hash table)
//                      module->selectedLineID = (int) lineDeets->getLineIndex();
//                      module->selectedPartID = (int) lineDeets->getPartIndex();
//                      module->hcOutput.setValue(0,0);
//                      module->hcPartnerLineID = 0;
//                      module->hcPartnerPartID = 0;
//                      //printf(" hcValues cleared LINE %d\n",__LINE__);
//                      module->show = true;
//                      module->integrate = true;
//                      //Update TopoDesign Scene render with new design information
//                      //theMsg->printf("Computing augmented connection:");
//                      if(module->computeAugmentedConnection() ) {
//                        module->topoDesignSelected = true;
//                        module->renderWWTopoDesignScene(false);
//                      }
//                  } else {
//                      //If no intersection, show design as incorrect, but still show dV line
//                      //Indicate that we don't have valid design point
//                      module->topoDesignSelected = false;
//                      //DO NOT clear the scene, just leave last design up
//                  }
//
//                  //Render the simple scene (just the dv line with points and no text)
//                  module->renderWWTempScene(false);
//
//              } //End Design start check
//              //If NOT started or button not down, then we actually follow the very first conditional (not wwDesign)
//            }
//    
//            //Call the compute function to run the rendering call and function evaluations
//            module->compute();
//        } //End Mouse intersects something
//      } //End Motion only check
//    } //End Location events
//
//    
//    //----------------------------------------------------------------------------------------
//    // Mouse Button Events - Only when no input is detected 
//    //----------------------------------------------------------------------------------------
//    /*theMsg->printf("Debug: hasInput = %d , Button1Press = %d , wwDesignStarted = %d ",
//                  hasInput, module->buttonPressed , wwDesignStarted);
//    theMsg->printf("Debug: button1Down = %d, EventIsLoc2Type = %d, EventIsButtonType = %d",
//                   module->buttonPressed, isLoc2Event, isButtonEvent);*/
//    //Registering a click indicates that we are initializing the TopoDesign selection between opposing manifold types
//    if (!hasInput && isButtonEvent) {
//        //Map to button events 
//        const SoMouseButtonEvent *mbEvent = (SoMouseButtonEvent*) eventCallback->getEvent();
//        //IF DESIGN IS NOT STARTED--------------------------------------------------------
//        if(!wwDesignStarted) {
//            //If pressed, lets see if we intersect the right node to start design
//            if (SO_MOUSE_PRESS_EVENT(mbEvent, BUTTON1) ) { //Clicking mouse button
//                //Debug
//                //theMsg->printf("%s : Button 1 Press - Looking for Target...",module->getName());
//                module->buttonPressed = true;
//                //Force Topology design to clear the rendering of a selected design
//                module->renderWWTopoDesignScene(true); //At new click
//                //Force that no design is selected
//                module->wwTopoDesignON = false; //Will change if valid point is selected
//                module->topoDesignSelected = false;
//                
//                //Only proceed if you intersect something relevant
//                const SoPickedPoint * pp = eventCallback->getPickedPoint();
//                if (pp!=NULL) {
//                  SbVec3f thePoint = pp->getPoint();
//                  //Look for the right manifold type 
//                  bool wuDesign = true; //Looking for unstable to start
//                  if(module->portDesignOpt.getValue() == 5) wuDesign=false; //Looking for stable to start
//                  //Check the scene to ensure intersecting the intended target
//                  const SoPickedPointList ppList = eventCallback->getAction()->getPickedPointList();
//                  bool containsWTarget = false; int wPathIdx = 0;
//                  //Loop through list to see if it holds the correct node
//                  for(int i=0;i<ppList.getLength();i++) {
//                    if(!containsWTarget) {
//                      if (wuDesign) {
//                        //Look for unstable manifold node
//                        containsWTarget = ppList[i]->getPath()->containsNode( module->wuTarget );
//                      } else {
//                        //Look for stable manifold node
//                        containsWTarget = ppList[i]->getPath()->containsNode( module->wsTarget );
//                      }
//                      if (containsWTarget) wPathIdx = i;
//                    }
//                  }
//                  //If target is selected:
//                  if(containsWTarget) {
//                    //Select the manifold point and store information
//                    SoLineDetail *lineDeets = (SoLineDetail*) ppList[wPathIdx]->getDetail();
//                    module->selectedUnstable = wuDesign;
//                    //Gather point information from line details (for hash table)
//                    module->selectedLineID = (int) lineDeets->getLineIndex();
//                    module->selectedPartID = (int) lineDeets->getPartIndex();
//                    module->pickedPoint = thePoint;
//                    //Disable all HC design options
//                    module->hcOutput.setValue(0,0);
//                    module->hcPartnerLineID = 0;
//                    module->hcPartnerPartID = 0;
//                    //Tell module to show and compute arcs if enabled
//                    module->show = true;
//                    module->integrate = true;
//                    //Initiate the WtoW Design (wwTopoDesignON = true)
//                    module->startAugmentedConnection();
//                    //Indicate that we have wtowDesign
//                    //module->topoDesignSelected = false; //Second point not defined
//                    //Render Temporary Topology Design Scene until button release
//                    module->renderWWTempScene(false);
//                    //theMsg->printf("%s : Target Found.  wwTopoDesignON = %d",module->getName(),module->wwTopoDesignON);
//                  }
//                  else {
//                    //theMsg->printf(" Not selected correct initial manifold (%s)",(wuDesign)?"Wu":"Ws");
//                  }
//                } //End pickedPoint intersects something
//            }
//            //And what if we release while wwDesign is not started 
//            if (SO_MOUSE_RELEASE_EVENT(mbEvent, BUTTON1)) {
//                //Clear the design as invalid 
//                module->renderWWTempScene(true); //Clear scene 
//                module->buttonPressed = false;
//                module->wwTopoDesignON = false;
//                //Basically does nothing...
//            }
//        }
//        //IF DESIGN IS STARTED--------------------------------------------------------
//        else {
//            //Release event finalizes the design
//            if (SO_MOUSE_RELEASE_EVENT(mbEvent, BUTTON1) ) {
//                //Clear the temporary selection scene
//                module->topoDesignSelected = false;
//                module->buttonPressed = false;
//                module->wwTopoDesignON = false;
//                module->renderWWTempScene(true);
//             
//                //Freezes a valid design 
//                //  Note: already computed during move events and data is assigned in
//                //  module->augHConn member.
//                
//                //Invalid design is not saved and everything is reset
//                if ( !(module->augHConn.isValid) ) {
//                  module->augHConn.reset();
//                  //Done running design
//                  module->wwTopoDesignON = false;
//                  //Clear the topology design scene
//                  module->renderWWTopoDesignScene(true);
//                  //theMsg->printf("%s: Selected Manifold to Manifold design is invalid! Clearing...",module->getName());
//                }
//
//            }
//            //If select again, this resets, but a release has to occur first.
//        }
//        //Render any changes in main scene or function evaluations
//        module->compute();
//    } //End button press check
//    
//
//}
//
///// Callback for button presses (keyboard)
//void DisplayManifolds_keyboardCB(void *userData, SoEventCallback *eventCallback)
//{
//    //Call this module
//    DisplayManifolds* module = (DisplayManifolds*) userData;
//    //Interface for QStrings
//    //QStringInterface qi;
//    
//    //The event
//    const SoEvent *event = eventCallback->getEvent();
//    //Run at keyboard button press
//    if (event->isOfType(SoKeyboardEvent::getClassTypeId()) ) {
//      //Output Trajectory 
//      if (event->wasCtrlDown() && SO_KEY_PRESS_EVENT(event,SPACE)) {
//        //Press button
//        module->simData.hit(2);
//        module->compute();
//        module->simData.clearHitState(2);
//      }
//      //Output HxLineSet
//      if (event->wasCtrlDown() && event->wasShiftDown() && SO_KEY_PRESS_EVENT(event,SPACE)) {
//        //Press button
//        module->simData.hit(0);
//        module->compute();
//        module->simData.clearHitState(0);
//      }
//      if (SO_KEY_RELEASE_EVENT(event,SPACE)) {
//        //Clear button press states
//        module->simData.clearHitState(0);
//        module->simData.clearHitState(2);
//      }
//      //Output Returns (and sub-iterates) to HxCluster object 
//      if (SO_KEY_PRESS_EVENT(event,R)) {
//        //Press button 
//        module->simData.hit(3);
//        module->compute();
//        module->simData.clearHitState(3);
//      }
//      if (SO_KEY_RELEASE_EVENT(event,R)) {
//          //Clear hit state
//          module->simData.clearHitState(3);
//      }
//      //Heteroclinic/Homoclinic Extraction if enabled
//      if (event->wasCtrlDown() && SO_KEY_PRESS_EVENT(event,H)) {
//          //Press button (which has logic embedded)
//          module->hcOutput.setValue(2,1); //Button press 
//          module->compute();
//          module->hcOutput.setValue(2,0); //Reset button
//      }
//      
//      //Run compute to execute button presses 
//      module->compute();
//    }
//}
//
////Integrate a manifold trajectory from selected (lineID,segID) pair 
//void DisplayManifolds::integrateManifoldArc(const bool useLinearParam)
//{
    //Check for input data
    ManifoldData *mData = (ManifoldData*) portData.source();
    if(!mData) {
      theMsg->printf("%s: Cannot propagate without valid manifold data",__FILE__);
      return;
    }
    
    //Check for Connected FixedPointData
    FixedPointData *fpData = (FixedPointData*) mData->fpxDataConnection.source();
    
    if(!fpData) {
      theMsg->printf("%s: Invalid FixedPointData connected to input!",__FILE__);
      return;
    }
    
    //Convert segment, line index, and parameter into manifold and state
    LineSegPair pickedManifold( selectedLineID, selectedPartID );
    LineSegPair hcPartnerManifold( hcPartnerLineID, hcPartnerPartID );
    //LineID,SegID -> ManifoldID
    IntPair manSegPair;
    IntPair hcManSegPair;
    if(selectedUnstable) {
      //Get (manifoldID,segmentID)
      manSegPair = uLinesetToManifold[pickedManifold];
    } else {
      manSegPair = sLinesetToManifold[pickedManifold];
    }
    //Heteroclinic piece
    if(hcOutput.getValue(0)==1) {
      if(selectedUnstable) {
        //Get (manifoldID,segmentID) for HC partner
        hcManSegPair = sLinesetToManifold[hcPartnerManifold];
      } else {
        hcManSegPair = uLinesetToManifold[hcPartnerManifold];
      }
    }
    
    //Convert the picked point into map space with the projector
    SbVec3f mapPoint; //(x,xdot,ydot)
    mapPoint = mapSpaceProjector.worldToMap(horzMenu,vertMenu,horzScale,vertScale,pickedPoint);
    nvis::vec2 xApprox(mapPoint[0],mapPoint[1]);
    //Get the Manifold and Segment 
    int mID = manSegPair.first; int segID = manSegPair.second;
    const ManifoldData::ManifoldType& theManifold = mData->theManifoldData.mapManifolds[mID];
    const ManifoldData::ManifoldSeg& theSeg = theManifold.segments[segID];
    
    //Compute the linear parameter for selected point
    if (!useLinearParam) { 
      //This is actually the default option.
      selectedLinearParameter = theSeg.getLinearParam(xApprox);
    } //Otherwise, use current value of selectedLinearParameter
    //Separate function for heteroclinic connection : Compute the exact intersection...
    
    
    //Look up the Progenitor State 
    //vec6 pState = theManifold.getPSOrbitState(segID,selectedLinearParameter);
    //progStatePoint.setValue( pState[0], pState[1], pState[2] );

    //Debug
    /*std::cout << "PickedManifold:  \n [lineIdx=" << selectedLineID << ",partIdx=" << selectedPartID
              << "] = [mID=" << manSegPair.first << " ("
              << ( (selectedUnstable)?"Unstable":"Stable" )
              << "), segID=" << manSegPair.second << "]\n";
    std::cout << " tau = " << selectedLinearParameter << "\n";
    //Inverse lookup
    std::cout << "Available Manifolds:\n";
    int numManifolds = (int) mData->theManifoldData.mapManifolds.size();
    for (int m=0;m<numManifolds;m++) {
      bool isUnstable = mData->theManifoldData.mapManifolds[m].isForward();
      std::cout << " [" << m << "] = " << ( (isUnstable)? "Unstable" : "Stable" ) << "\n";
    }*/
    
    //Debugging:
    selectedParts.setValue(0,selectedLineID);
    selectedParts.setValue(1,selectedPartID);
    selectedManSeg.setValue(0,mID);
    selectedManSeg.setValue(1,segID);
    hcPartnerParts.setValue(0,hcPartnerLineID);
    hcPartnerParts.setValue(1,hcPartnerPartID);
    if(hcOutput.getValue(0) == 1) {
      hcPartnerManSeg.setValue(0,hcManSegPair.first);
      hcPartnerManSeg.setValue(1,hcManSegPair.second);
    } else {
      hcPartnerManSeg.setValue(0,-1);
      hcPartnerManSeg.setValue(1,-1);
    }
    QString tauString;
    tauString.setNum( selectedLinearParameter );
    selectedParamText.setValue( tauString.toLocal8Bit().data() );
    
    
    //Utilizing a linear parameter is more accurate for manifold seeding 
    //than just a point (float) picked from manifold
    //    ManifoldID,segmentID,linParam -> map state (x,xdot) 
    nvis::vec2 x0 = theSeg.getPoint(selectedLinearParameter);
    
    //Evaluate length of propagation (based on depth)
    int period = theManifold.getPeriod();
    //Forward or backward determined by upstream/downstream and manifold type
    int maxDepth = 0, currentDepth = 0;
    try {
        maxDepth = theManifold.getDepth();
        currentDepth = theManifold.getDepth(segID);
    } catch (const std::exception& e) {
        theMsg->printf("%s: Error while searching for depth %s ",__FILE__,e.what() );
        return;
    }
    selectedPropInfo.setValue(0,currentDepth);
    //Upstream and downstream periods determined by depth
    int pDownstream = period * (maxDepth - currentDepth);
    //Note: To get to depth = 0
    int pUpstream = -period * (currentDepth);
    
    //Switch between propagation types:
    //AND Gather the iterates utilizing data lookup from the manifold
    int p = pUpstream;
    bool upstreamSim = false;
    QString propTypeStr("");
    uIterates.clear(); uItTimes.clear();
    dIterates.clear(); dItTimes.clear();
    switch (simDirection.getValue()) {
       case 0 :
        // 1) Propagate the upstream portion
        propTypeStr = "Upstream";
        p = pUpstream;
        upstreamSim = true;
        mData->theManifoldData.mapManifolds[mID].getUpstreamPoints(segID,selectedLinearParameter,uIterates);
        //Data Display
        if(simOptions.getValue(2)==1) {
          //May need times 
          theMsg->printf("Upstream points on selected manifold arc:");
          for(int i=0;i<(int)uIterates.size();i++) {
            theMsg->printf(" %d : [ %g %g ] ",i, uIterates[i][0], uIterates[i][1] );
          }
        }
        break;
       case 1:
        // 2) Propagate the downstream portion
        propTypeStr = "Downstream";
        p = pDownstream;
        mData->theManifoldData.mapManifolds[mID].getDownstreamPoints(segID,selectedLinearParameter,dIterates);
        //Data display
        if (simOptions.getValue(2)==1) {
          //May need times
          theMsg->printf("Downstream points on selected manifold arc:");
          for(int i=0;i<(int)dIterates.size();i++) {
            theMsg->printf(" %d : [ %g %g ] ",i, dIterates[i][0], dIterates[i][1] );
          }
        }
        break;
       case 2:
        // 3) Propagate both upstream and downstream
        propTypeStr = "Both";
        p = pUpstream; //Show upstream
        upstreamSim = true;
        //Note, will do both up and downstream propagation
        mData->theManifoldData.mapManifolds[mID].getUpstreamPoints(segID,selectedLinearParameter,uIterates);
        mData->theManifoldData.mapManifolds[mID].getDownstreamPoints(segID,selectedLinearParameter,dIterates);
        //May need times
        break;
    }
    selectedPropInfo.setLabel(propTypeStr);
    selectedPropInfo.setValue(1,p);
    
    //Clear the current data for actual trajectories
    uStates.clear(); uTimes.clear();
    dStates.clear(); dTimes.clear();
    
    //Also display the selected Orbit and fixed point indexes
    if(simOptions.getValue(2)) {
      theMsg->printf(" Selected From: OrbitID = %d  FixedPoint = %d", 
                     mData->theManifoldData.mapManifolds[mID].fpdOrbitIdx,
                     mData->theManifoldData.mapManifolds[mID].fpdPointIdx);
    }
    
    //Stop here if we don't need to actually propagate
    if (simOptions.getValue(0)==0 && returnDisplayOpts.getValue(2)==0) return;
    
    
    //Build integration objects
    rhs_type* theEOMs = new rhs_type(C,mup);  //The right-hand side:  CR3BP EoMs
    section_type *theSection = new section_type(*theEOMs); //Hyperplane: y=0 by planar_section
    map_type pMap(*theEOMs, *theSection); //Map object - works like an engine
    pMap.setPrecision( 1.e-12 ); //integration tolerance
    //Temporary objects for mapping
    std::vector<double> psTimes;
    std::vector<vec6> psStates;
    
    //Always Propagate/Store UPSTREAM first
    if (upstreamSim) {
        //Upstream simulations always propagate to Progenitor State
        std::vector<vec6>::iterator it;
        
        //Run the propagation for upstream points
        if(returnDisplayOpts.getValue(2) == 1) {
            //We have to use the sub-return version of the propagation
            mData->theManifoldData.mapManifolds[mID].propagateUpstream(
               pMap,segID,selectedLinearParameter, psStates, psTimes,
               uIterates, uItTimes, uSubIts, uSubItsT);
        } else {
            //Employ the version that just returns states (as vec6's)
            mData->theManifoldData.mapManifolds[mID].propagateUpstream(
                pMap,segID,selectedLinearParameter, psStates, psTimes);
        }
        if ((int)psStates.size() < 1) {
            theMsg->printf("%s: Map error - %s", __FILE__, "Failed UPSTREAM propagation");
            return;
        }
        
        //Store states to a State
        //vec42 y = theSection->unproject(x0);
        //uStates.push_back(State(y[0],y[1],y[2],y[3],y[4],y[5]));
        //uTimes.push_back(0.0);
        for(int i=0; i<(int)psStates.size();i++) {
            vec6& yps = psStates[i];
            uStates.push_back( State(yps[0],yps[1],yps[2],yps[3],yps[4],yps[5]) );
            uTimes.push_back( psTimes[i] );
        }
        //theMsg->printf("Integration called:  Num integration states = %d", out_states.size());

        //Display state and final point
        vec42 icState = theSection->unproject(x0);
        if(simOptions.getValue(2)) {
            theMsg->printf(" State of Selection: [%g %g %g %g %g %g]",
                            icState[0],icState[1],icState[2],
                            icState[3],icState[4],icState[5]);
            theMsg->printf(" Propagation: Upstream");
            theMsg->printf(" Depth trace:  MaxDepth = %d Current Depth = %d", maxDepth, currentDepth);
            theMsg->printf(" Arc Prop Returns (p) = %d with total time %g", pUpstream, uTimes.back());
        
        }
    }
    
    //Clear data (not to overwrite)
    psStates.clear(); psTimes.clear();
    
    // Propagate/Store the DOWNSTREAM information
    if (simDirection.getValue() == 0 || simDirection.getValue() == 2) {
        //Call propagation from manifold component:
        if(returnDisplayOpts.getValue(2) == 1) {
            //Utilize the sub-return version of the propagation :
            mData->theManifoldData.mapManifolds[mID].propagateDownstream(
                pMap, segID, selectedLinearParameter, psStates, psTimes,
                dIterates, dItTimes, dSubIts, dSubItsT);
        } else {
            //Employ just the internal states version (reinterpolate between manifold crossings on map)
            mData->theManifoldData.mapManifolds[mID].propagateDownstream(
                pMap, segID, selectedLinearParameter, psStates, psTimes);
        }
        if((int)psStates.size() <1 ) {
            theMsg->printf("%s: Map exception - %s", __FILE__, "Failed DOWNSTREAM propagation");
            return;
        }
    
        //Store states to a State
        vec42 y = theSection->unproject(x0);
        //dStates.push_back(State(y[0],y[1],y[2],y[3],y[4],y[5]));
        //dTimes.push_back(0.0);
        for(int i=0;i<(int)psStates.size();i++) {
          vec6& yps = psStates[i];
          dStates.push_back(State(yps[0],y[1],yps[2],yps[3],yps[4],yps[5]));
          dTimes.push_back(psTimes[i]);
        }
        
        //Display state and final point
        vec42 icState = theSection->unproject(x0);
        if(simOptions.getValue(2)) {
            theMsg->printf(" State of Selection: [%g %g %g %g %g %g]",
                            icState[0],icState[1],icState[2],
                            icState[3],icState[4],icState[5]);
            theMsg->printf(" Propagation: Downstream");
            theMsg->printf(" Depth trace:  MaxDepth = %d Current Depth = %d", maxDepth, currentDepth);
            theMsg->printf(" Arc Prop Returns (p) = %d with total time %g", pDownstream, dTimes.back());
        
            //std::vector<return_type>::reverse_iterator rit;
            //rit = rTypeIters.rbegin();
            //nvis::vec2 mapDisp = rit->x - x0;
            //theMsg->printf(" Map Displacement = P^p(x) - x = [%g %g]  Mag=%g",
            //               mapDisp[0],mapDisp[1],nvis::norm(mapDisp));
        }
    }
    

}


//Initialization
// HX_INIT_CLASS(DisplayManifolds,HxModule)

//==========================================================================================
//Constructor::
// DisplayManifolds::DisplayManifolds() :
//     HxModule(ManifoldData::getClassTypeId()),
//     info(this,"Info",QApplication::translate("DisplayManifolds","Info")),
//     initStateConnection(this,"InputState",QApplication::translate("DisplayManifolds","InputState"),HxInitialState::getClassTypeId() ),
//     portTabs(this,"Category",QApplication::translate("DisplayManifolds","Category"),5),
//     system(this,"System",QApplication::translate("DisplayManifolds","System"),5),
//     gravParam(this,"Mu",QApplication::translate("DisplayManifolds","Mu")),
//     portJacobiConstant(this,"Jacobi",QApplication::translate("DisplayManifolds","Jacobi")),
//     portXbounds(this,"Bounds(X)",QApplication::translate("DisplayManifolds","Bounds(X)"),2),
//     portYbounds(this,"Bounds(Y)",QApplication::translate("DisplayManifolds","Bounds(Y)"),2),
//     sep0(this,"Sep",QApplication::translate("DisplayManifolds","Sep")),
//     displayOptions(this,"Options",QApplication::translate("DisplayManifolds","Options"),4),
//     horzMenu(this,"HorzAxis",QApplication::translate("DisplayManifolds","HorzAxis"),6),
//     vertMenu(this,"VertAxis",QApplication::translate("DisplayManifolds","VertAxis"),6),
//     horzScale(this,"HorzScale",QApplication::translate("DisplayManifolds","HorzScale")),
//     vertScale(this,"VerScale",QApplication::translate("DisplayManifolds","VertScale")),
//     eigenLength(this,"EigLength",QApplication::translate("DisplayManifolds","EigLength")),
//
//     sep1(this,"Sep",QApplication::translate("DisplayManifolds","Sep")),
//     dispInfo(this,"Info",QApplication::translate("DisplayManifolds","Info")),
//     colors(this,"Colors",QApplication::translate("DisplayManifolds","Colors"),3),
//     stableColormap(this,"Stable"),
//     unstableColormap(this,"Unstable"),
//     portLineWidth(this,"Line Width",QApplication::translate("DisplayManifolds","Line Width")),
//     portPointSize(this,"Point Size",QApplication::translate("DisplayManifolds","Point Size")),
//     portMarkerScale(this,"Marker Scale",QApplication::translate("DisplayManifolds","Marker Scale")),
//     sep2(this,"Sep",QApplication::translate("DisplayManifolds","Sep")),
//     filterInfo(this,"Filter",QApplication::translate("DisplayManifolds","Filter")),
//     filterOptions(this,"FOpts",QApplication::translate("DisplayManifolds","FOpts"),2),
//     filterType(this,"FType",QApplication::translate("DisplayManifolds","FType"),3),
//     filterOrbitID(this,"FOrbitID",QApplication::translate("DisplayManifolds","FOrbitID")),
//     filterFixedPoint(this,"FFPID",QApplication::translate("DisplayManifolds","FFPID")),
//
//     selectionInfo(this,"Selection",QApplication::translate("DisplayManifolds","Selection")),
//     simOptions(this,"SimOpt",QApplication::translate("DisplayManifolds","SimOpt"),3),
//     simDirection(this,"SimDir",QApplication::translate("DisplayManifolds","SimDir"),3),
//     simRevColor(this,"RCol",QApplication::translate("DisplayManifolds","RCol"),2),
//     simData(this,"SimData",QApplication::translate("DisplayManifolds","SimData"),4),
//     hcOutput(this,"HConn",QApplication::translate("DisplayManifolds","HConn")),
//     selectedParts(this,"Parts",QApplication::translate("DisplayManifolds","Parts"),2),
//     selectedManSeg(this,"ManSeg",QApplication::translate("DisplayManifolds","ManSeg"),2),
//     selectedPropInfo(this,"PropInfo",QApplication::translate("DisplayManifolds","PropInfo"),2),
//     selectedParamText(this,"Param",QApplication::translate("DisplayManifolds","Param")),
//     hcPartnerParts(this,"HCParts",QApplication::translate("DisplayManifolds","HCParts"),2),
//     hcPartnerManSeg(this,"HCManSeg",QApplication::translate("DisplayManifolds","HCManSeg"),2),
//     sep3(this,"Sep",QApplication::translate("DisplayManifolds","Sep")),
//     returnInfo(this,"RInfo",QApplication::translate("DisplayManifolds","RInfo")),
//     returnDisplayOpts(this,"ROpts",QApplication::translate("DisplayManifolds","ROpts"),3),
//     returnMarkerScale(this,"RScale",QApplication::translate("DisplayManifolds","RScale")),
//     trajOutInfo(this,"TrajHK",QApplication::translate("DisplayManifolds","TrajHK")),
//     trajOutInfo2(this,"LineHK",QApplication::translate("DisplayManifolds","LineHK")),
//     returnsOutInfo(this,"Returns",QApplication::translate("DisplayManifolds","Returns")),
//     hcOutInfo(this,"HCTrajs",QApplication::translate("DisplayManifolds","HCTrajs")),
//     sep4(this,"Sep",QApplication::translate("DisplayManifolds","Sep")),
//
//     funcInfo(this,"Funcs",QApplication::translate("DisplayManifolds","Funcs")),
//     funcOrbitID(this,"DataOrbID",QApplication::translate("DisplayManifolds","DataOrbID")),
//     funcManifoldID(this,"DataManID",QApplication::translate("DisplayManifolds","DataManID"),1),
//     clearFuncs(this,"Clear",QApplication::translate("DisplayManifolds","Clear"),2),
//     portOrbitIDList(this,"IDList",QApplication::translate("DisplayManifolds","IDList")),
//     extractFuncs(this,"Extract",QApplication::translate("DisplayManifolds","Extract"),4),
//
//     portDesignInfo(this,"Design",QApplication::translate("DisplayManifolds","Design")),
//     jcIS(this,"Jacobi0",QApplication::translate("DisplayManifolds","Jacobi0")),
//     portDesignOpt(this,"DesignOpt",QApplication::translate("DisplayManifolds","DesignOpt"),6),
//     portDesignDisplayOpts(this,"DesignDisp",QApplication::translate("DisplayManifolds","DesignDisp"),3),
//     portDeltaVUnits(this,"DesignUnits",QApplication::translate("DisplayManifolds","DesignUnits"),3),
//     designColors(this,"DColors",QApplication::translate("DisplayManifolds","DColors"),4),
//     portDesignFont(this,"DFont",QApplication::translate("DisplayManifolds","DFont")),
//     portDesignOutput(this,"DOut",QApplication::translate("DisplayManifolds","DOut"),2),
//     portWtoWDesignOutput(this,"WWDesign",QApplication::translate("DisplayManifolds","WWDesign"),3)
//
// {
//     //Default values of ports
//     info.setValue("Saddle manifold display (For y=0 planar map only!)");
//
//     //Tabs
//     portTabs.setLabel(0,"Problem");
//     portTabs.portAdd(0,&system);
//     portTabs.portAdd(0,&gravParam);
//     portTabs.portAdd(0,&portJacobiConstant);
//     portTabs.portAdd(0,&portXbounds);
//     portTabs.portAdd(0,&portYbounds);
//     portTabs.portAdd(0,&sep0);
//     portTabs.portAdd(0,&horzMenu);
//     portTabs.portAdd(0,&vertMenu);
//     portTabs.portAdd(0,&horzScale);
//     portTabs.portAdd(0,&vertScale);
//
//     //Display
//     portTabs.setLabel(1,"Display");
//     portTabs.portAdd(1,&displayOptions);
//     portTabs.portAdd(1,&eigenLength);
//     portTabs.portAdd(1,&sep1);
//     portTabs.portAdd(1,&dispInfo);
//     portTabs.portAdd(1,&colors);
//     portTabs.portAdd(1,&stableColormap);
//     portTabs.portAdd(1,&unstableColormap);
//     portTabs.portAdd(1,&portLineWidth);
//     portTabs.portAdd(1,&portPointSize);
//     portTabs.portAdd(1,&portMarkerScale);
//     portTabs.portAdd(1,&sep2);
//     portTabs.portAdd(1,&filterInfo);
//     portTabs.portAdd(1,&filterOptions);
//     portTabs.portAdd(1,&filterType);
//     portTabs.portAdd(1,&filterOrbitID);
//     portTabs.portAdd(1,&filterFixedPoint);
//     portTabs.portAdd(1,&sep3);
//     portTabs.portAdd(1,&returnInfo);
//     portTabs.portAdd(1,&returnDisplayOpts);
//     portTabs.portAdd(1,&returnMarkerScale);
//
//     //Selection
//     portTabs.setLabel(2,"Select");
//     portTabs.portAdd(2,&selectionInfo);
//     portTabs.portAdd(2,&simOptions);
//     portTabs.portAdd(2,&simDirection);
//     portTabs.portAdd(2,&simRevColor);
//     portTabs.portAdd(2,&simData);
//     portTabs.portAdd(2,&hcOutput);
//     portTabs.portAdd(2,&selectedParts);
//     portTabs.portAdd(2,&selectedManSeg);
//     portTabs.portAdd(2,&selectedPropInfo);
//     portTabs.portAdd(2,&selectedParamText);
//     portTabs.portAdd(2,&hcPartnerParts);
//     portTabs.portAdd(2,&hcPartnerManSeg);
//
//     //Functions
//     portTabs.setLabel(3,"Functions");
//     portTabs.portAdd(3,&funcInfo);
//     portTabs.portAdd(3,&selectedManSeg);
//     portTabs.portAdd(3,&funcOrbitID);
//     portTabs.portAdd(3,&funcManifoldID);
//     portTabs.portAdd(3,&clearFuncs);
//     portTabs.portAdd(3,&portOrbitIDList);
//     portTabs.portAdd(3,&extractFuncs);
//     portTabs.portAdd(3,&sep2);
//     portTabs.portAdd(3,&filterInfo);
//     portTabs.portAdd(3,&filterOptions);
//     portTabs.portAdd(3,&filterType);
//     portTabs.portAdd(3,&filterOrbitID);
//     portTabs.portAdd(3,&filterFixedPoint);
//     portTabs.portAdd(3,&sep3);
//     portTabs.portAdd(3,&trajOutInfo);
//     portTabs.portAdd(3,&trajOutInfo2);
//     portTabs.portAdd(3,&returnsOutInfo);
//     portTabs.portAdd(3,&hcOutInfo);
//     portTabs.portAdd(3,&sep4);
//
//     //Topology Design Module
//     portTabs.setLabel(4,"TopoDesign");
//     portTabs.portAdd(4,&portDesignInfo);
//     portTabs.portAdd(4,&jcIS);
//     portTabs.portAdd(4,&portDesignOpt);
//     portTabs.portAdd(4,&portDesignDisplayOpts);
//     portTabs.portAdd(4,&portDeltaVUnits);
//     portTabs.portAdd(4,&designColors);
//     portTabs.portAdd(4,&portDesignFont);
//     portTabs.portAdd(4,&portDesignOutput);
//     portTabs.portAdd(4,&portWtoWDesignOutput);
//
//
//
//     //System port
//     sys = new CR3BPsystem();
//     sys->buildMenus(system,gravParam);
//     mup = sys->mup;
//     portJacobiConstant.setValue("2.96");
//     C = 2.96;
//     portXbounds.setLabel(0,"Min");
//     portXbounds.setLabel(1,"Max");
//     portXbounds.setValue(0,-100);
//     portXbounds.setValue(1,100);
//     portYbounds.setLabel(0,"Min");
//     portYbounds.setLabel(1,"Max");
//     portYbounds.setValue(0,-2.55);
//     portYbounds.setValue(1,2.55);
//
//
//
//
//     //Mapping to world space (likely projection to Rotating frame)
//     mapSpaceProjector.buildMenu(horzMenu);
//     mapSpaceProjector.buildMenu(vertMenu);
//     vertMenu.setValue(2);
//     horzScale.setValue(1);
//     vertScale.setValue(0.16); //Standard scale in EM image
//
//
//     //Eigen Space (Vectors)
//     displayOptions.setLabel(0,"EigVec");
//     displayOptions.setLabel(1,"FPX");
//     displayOptions.setLabel(2,"MapDiscont");
//     displayOptions.setLabel(3,"Manifolds");
//     displayOptions.setValue(0,0);
//     displayOptions.setValue(1,0);
//     displayOptions.setValue(2,0);
//     displayOptions.setValue(3,1);
//     eigenLength.setMinMax(0.01,1);//NonDim Length
//     eigenLength.setValue(0.2);
//
//     //--- Display ---
//     dispInfo.setValue("Display Settings");
//     McColor *green = new McColor;
//     green->setValue(0.5,1.0,0.5);
//     McColor *khaki = new McColor;
//     khaki->setValue(195./255.,176./255.,145./255.);
//     McColor *red = new McColor;
//     red->setValue(1.0,0.2,0.2);
//     McColor *smalt = new McColor;
//     smalt->setValue(0,51./255.,153./255.);
//     colors.setLabel(0,"FixPts");
//     colors.setColor(*green,0);
//     colors.setLabel(1,"Stable");
//     colors.setColor(*smalt,1);
//     colors.setLabel(2,"Unstable");
//     colors.setColor(*red,2);
//     stableColormap.setLabel("Stable");
//     stableColormap.setDefaultColor( SbColor(smalt->r,smalt->g,smalt->b) );
//     unstableColormap.setLabel("Unstable");
//     unstableColormap.setDefaultColor( SbColor(red->r,red->g,red->b) );
//     portLineWidth.setMinMax(0,10);  //Only goes from 0 to 10
//     portLineWidth.setValue(2);
//     portPointSize.setMinMax(0,10);
//     portPointSize.setValue(2);
//     portMarkerScale.setMinMax(0.01,10);
//     portMarkerScale.setValue(1);
//
//
//     //Simple display filters
//     filterInfo.setValue("Simple filters for displaying components.");
//     filterOptions.setLabel(0,"SingleOrbit");
//     filterOptions.setLabel(1,"SingleFP");
//     filterType.setLabel(0,"Stable");
//     filterType.setLabel(1,"Unstable");
//     filterType.setLabel(2,"Both");
//     filterType.setValue(2);
//     filterOrbitID.setMinMax(0,0);
//     filterOrbitID.hide();
//     filterFixedPoint.setMinMax(0,0);
//     filterFixedPoint.hide();
//
//     //--- Select ----
//     selectionInfo.setValue("-----");
//     simOptions.setLabel(0,"Arc");
//     simOptions.setLabel(1,"Returns");
//     simOptions.setLabel(2,"Data");
//     //simOptions.setLabel(3,"PState");
//     simOptions.setValue(1,1);
//     simDirection.setLabel(0,"Upstream");
//     simDirection.setLabel(1,"Downstream");
//     simDirection.setLabel(2,"Both");
//     simRevColor.setLabel(0,"Same");
//     simRevColor.setLabel(1,"Compliment");
//     simData.setLabel(0,"Lineset");
//     simData.setLabel(1,"Data");
//     simData.setLabel(2,"Traj");
//     simData.setLabel(3,"Returns");
//     hcOutput.insertCheckBox(0,"Found");
//     hcOutput.insertColorButton(1,"");
//     hcOutput.setColor(1,(*green));
//     hcOutput.insertPushButton(2,"OutputTraj");
//     hcOutput.insertPushButton(3,"OutputData");
//
//     //--- Debug Selection ---
//     selectedParts.setLabel(0,"LineID");
//     selectedParts.setLabel(1,"Part/SegID");
//     selectedParts.setValue(0,-1);
//     selectedParts.setValue(1,-1);
//     selectedManSeg.setLabel(0,"ManifoldID");
//     selectedManSeg.setLabel(1,"SegmentID");
//     selectedManSeg.setValue(0,-1);
//     selectedManSeg.setValue(1,-1);
//     selectedPropInfo.setLabel(0,"Depth");
//     selectedPropInfo.setLabel(1,"Iterates");
//     selectedPropInfo.setValue(0,0);
//     selectedPropInfo.setValue(1,0);
//     selectedParamText.setValue("0");
//     hcPartnerParts.setLabel(0,"LineID");
//     hcPartnerParts.setLabel(1,"Part/SegID");
//     hcPartnerParts.setValue(0,-1);
//     hcPartnerParts.setValue(1,-1);
//     hcPartnerManSeg.setLabel(0,"ManifoldID");
//     hcPartnerManSeg.setLabel(1,"SegmentID");
//     hcPartnerManSeg.setValue(0,-1);
//     hcPartnerManSeg.setValue(1,-1);
//
//     //Return display options
//     returnInfo.setValue("Display options for selected arc returns");
//     returnDisplayOpts.setLabel(0,"Selection");
//     returnDisplayOpts.setLabel(1,"Returns");
//     returnDisplayOpts.setLabel(2,"SubReturns");
//     returnDisplayOpts.setValue(0,1);
//     returnDisplayOpts.setValue(1,1);
//     returnMarkerScale.setMinMax(0.01,5);
//     returnMarkerScale.setValue(1);
//
//     //Hot-Key information
//     trajOutInfo.setValue("Ctrl+Spacebar : Output Trajectory from selection");
//     trajOutInfo2.setValue("Ctrl+Shift+Spacebar : Output HxLineSet from selection");
//     returnsOutInfo.setValue("R : Output HxCluster for returns (and sub-returns if marked)");
//     hcOutInfo.setValue("Ctrl+H : Output Trajectory objects for heteroclinic/homoclinic selection");
//
//     //Functions|Ops for manipulating input data and new objects
//     funcInfo.setValue("Data manipulation functions: Clearing and Extracting");
//     funcOrbitID.setMinMax(0,0);
//     funcManifoldID.setValue(0);
//     clearFuncs.setLabel(0,"Orbit");
//     clearFuncs.setLabel(1,"Manifold");
//     portOrbitIDList.setValue("");
//     extractFuncs.setLabel(0,"Single");
//     extractFuncs.setLabel(1,"AddID");
//     extractFuncs.setLabel(2,"RemoveID");
//     extractFuncs.setLabel(3,"FromList");
//
//
//     //Topology Design Options
//     portDesignInfo.setValue("Map Topology Design Interface");
//     jcIS.setValue("2.96");
//     portDesignOpt.setLabel(0,0,"State to Ws (F)");
//     portDesignOpt.setLabel(0,1,"Wu to State (F)");
//     portDesignOpt.setLabel(0,2,"State to Wu (B)");
//     portDesignOpt.setLabel(0,3,"Ws to State (B)");
//     portDesignOpt.setLabel(0,4,"Wu to Ws (F)");
//     portDesignOpt.setLabel(0,5,"Ws to Wu (B)");
//     portDesignOpt.setValue(0,4); //Assume no input.
//     portDesignDisplayOpts.setLabel(0,"dVText");
//     portDesignDisplayOpts.setValue(0,1);
//     portDesignDisplayOpts.setLabel(1,"OrbitID");
//     portDesignDisplayOpts.setValue(1,1);
//     portDesignDisplayOpts.setLabel(2,"dxText");
//     portDesignDisplayOpts.setValue(2,1);
//     portDeltaVUnits.setLabel(0,"ND");
//     portDeltaVUnits.setLabel(1,"km/s");
//     portDeltaVUnits.setLabel(2,"m/s");
//     portDeltaVUnits.setValue(2);
//     McColor *ygreen = new McColor;
//     ygreen->setValue(173./255., 255./255., 47./255.);
//     McColor *darkOrchid = new McColor;
//     darkOrchid->setValue(154./255., 50./255., 205./255.);
//     designColors.setLabel(0,"InitialState");
//     designColors.setLabel(1,"deltaV");
//     designColors.setLabel(2,"dVLine");
//     designColors.setLabel(3,"OrbitID");
//     designColors.setColor(*darkOrchid,0);
//     designColors.setColor(*ygreen,1);
//     designColors.setColor(McColor(0.0,0.0,0.0),2);
//     McColor  *ired = new McColor;
//     ired->setValue(238./255., 99./255., 99./255.);
//     designColors.setColor(*ired,3);
//     portDesignFont.setColorButtonVisible(true);
//     portDesignFont.setFontColor(McColor(0.0,0.0,0.0));
//     portDesignOutput.setLabel(0,"EndState");  //For the next design node...
//     QString dvButtonStr = QString(QChar(0x0394)) + "V";
//     portDesignOutput.setLabel(1,dvButtonStr);
//     portWtoWDesignOutput.setLabel(0,"Objects");
//     portWtoWDesignOutput.setLabel(1,"Print");
//     portWtoWDesignOutput.setLabel(2,"Clear");
//     deltaV.setValue(0.0,0.0,0.0);
//     hasInput = false;
//     buttonPressed = false;
//     topoDesignSelected = false;
//     wwTopoDesignON = false;
//
//
//     //Constructor for Scene graph (new SoNode)
//     scene = new SoSeparator;
//     trajScene = new SoSeparator;
//     tdScene = new SoSeparator;
//     wwdvScene = new SoSeparator;
//
//     //Event callback for interaction with this object (propagate arc)
//     eventCB = new SoEventCallback;
//     eventCB->addEventCallback(SoMouseButtonEvent::getClassTypeId(),DisplayManifolds_clickMotionCB,this);
//     eventCB->addEventCallback(SoLocation2Event::getClassTypeId(),DisplayManifolds_clickMotionCB,this);
//     eventCB->addEventCallback(SoKeyboardEvent::getClassTypeId(),DisplayManifolds_keyboardCB,this);
//
//     //Need eventCB at root-level because we have to grab from multiple scenes
//     root  = (SoSeparator*) theController->getSceneGraph(0);
//     initCB = false;
//
//     pickedPoint.setValue(0.0,0.0,0.0);
//     progStatePoint.setValue(0.0,0.0,0.0);
//     show = integrate = false;
//     selectedLinearParameter = 0.0;
//     selectedLineID = 0;
//     selectedPartID = 0;
//     hcPartnerLineID = 0;
//     hcPartnerPartID = 0;
//
//     //Render the manifolds to define the data hash tables
//     //renderManifolds();
//
// }


// //Destructor::
// DisplayManifolds::~DisplayManifolds()
// {
//     if (initCB) root->removeChild(eventCB);
//     hideGeom(wwdvScene);
//     wwdvScene->removeAllChildren();
//     hideGeom(tdScene);
//     tdScene->removeAllChildren();
//     hideGeom(scene);
//     scene->removeAllChildren();
//     hideGeom(trajScene);
//     trajScene->removeAllChildren();
// }


// //Compute:: -> Called when a port is modified
// void DisplayManifolds::compute()
// {
//     //Picking Interactions at root level
//     if (!initCB) {
//       root->insertChild(eventCB,0);
//       initCB=true;
//     }
//
//     //Check for input data
//     ManifoldData *mData = (ManifoldData*) portData.source();
//     if(!mData) {
//       theMsg->printf("%s: Invalid input data",__FILE__);
//       //Free up Jacobi port
//       portJacobiConstant.setSensitivity(true);
//       return;
//     }
//
//
//     //Check for Connected FixedPointData
//     FixedPointData *fpData = (FixedPointData*) mData->fpxDataConnection.source();
//     if(!fpData) {
//       theMsg->printf("%s: Invalid FixedPointData connected to input!",__FILE__);
//       return;
//     }
//
//
//     //Update Jacobi constant value if data is new
//     if(portData.isNew()) {
//       double setJC = fpData->getJacobiConstant();
//       QString strJC;
//       strJC.setNum(setJC,'g',14);
//       portJacobiConstant.setSensitivity(true);
//       portJacobiConstant.setValue(strJC.toLocal8Bit().data());
//       portJacobiConstant.setSensitivity(false);
//     }
//     // Get the input parameters from the user interface:
//     //mup = sys->compute(system,gravParam);
//     sys->computeWithJacobi(system,gravParam,portJacobiConstant,mup,C);
//
//     QStringInterface qInterface;
//
//     //Check for an Initial State object in input (restricts selection)
//     HxInitialState *iState = static_cast<HxInitialState*>(initStateConnection.source()); // Go STATE!!
//     if(!iState) { hasInput = false; } else hasInput = true;
//
//     //theMsg->printf("%s : InitialState object is %s",getName(),((hasInput)?"CONNECTED":"NULL"));
//     if(hasInput) {
//         if(initStateConnection.isNew()) {
//             //Set the Jacobi constant value
//             QString jcISStr;
//             jcISStr.setNum( sys->getJacobiConstant(iState->data.ic), 'g', 14);
//             jcIS.setSensitivity(true);
//             jcIS.setValue( jcISStr );
//             jcIS.setSensitivity(false);
//             //Set the design port (is this a good idea?)
//             portDesignOpt.setValue(0,0);
//             //Set the pickedPoint to state location (hopefully only the first time)
//             SbVec3f mapPoint(iState->data.ic.x,iState->data.ic.xd,0.0);
//             pickedPoint = mapSpaceProjector.mapToWorld(horzMenu,vertMenu,horzScale,vertScale,
//                                                        mapPoint);
//             //Have to render the topology design scene
//             renderTopoDesignScene();
//         }
//     }
//
//     //-----------------------------------------------------------------------------------
//     //Initialization
//     //-----------------------------------------------------------------------------------
//     //Data Structure Variables
//     nvis::bbox2 _bounds;
//     xavier::metric<double,2>  euclidean_metric;
//
//     //Set the min/max of each coord on plane
//     _bounds = nvis::bbox2(nvis::vec2(portXbounds.getValue(0), portYbounds.getValue(0)),
//                 nvis::vec2(portXbounds.getValue(1), portYbounds.getValue(1)));
//     //Note, need infinite bounds to not stop at edge of analysis boxes
//     //const double LARGE = std::numeric_limits<double>::max();
//     //_bounds.min() = nvis::vec2(-LARGE,-LARGE);
//     //_bounds.max() = nvis::vec2(LARGE,LARGE);
//
//     //Map space (i.e. plane) is not periodic in either direction
//     euclidean_metric.periodic()[0] = false;
//     euclidean_metric.periodic()[1] = false;
//     //Metric bounds?
//     //euclidean_metric.bounds() = _bounds;
//
//     //Adapt ports to new input data
//     if (portData.isNew() ) {
//         //Adjust the the filter IDs :
//         int numOrbits = fpData->getNumOrbits();
//         filterOrbitID.setMinMax(0,numOrbits-1);
//         funcOrbitID.setMinMax(0,numOrbits-1);
//     }
//     //Show/hide orbit filter ports
//     if (filterOptions.isNew()) {
//         if(filterOptions.getValue(0)==1 && filterOptions.getValue(1)==1) {
//             filterOrbitID.show();
//             filterFixedPoint.show();
//         } else if(filterOptions.getValue(0)==0 && filterOptions.getValue(1)==1) {
//             //Don't allow a particular fixed point to be selected
//             filterOptions.setValue(1,0);
//             filterOrbitID.hide();
//             filterFixedPoint.hide();
//         } else if(filterOptions.getValue(0)==1) {
//             filterOrbitID.show();
//         } else {
//             filterOrbitID.hide();
//             filterFixedPoint.hide();
//         }
//     }
//     //Modify number of orbits if new filtering orbit selected
//     if (filterOrbitID.isNew()) {
//         int orbitID = filterOrbitID.getValue();
//         const xavier::fixpoint& fp = fpData->getFixedPoint(orbitID,0);
//         int p = fp.K;
//         filterFixedPoint.setMinMax(0,p-1);
//     }
//
//     //Change options if C or input data is modified?
//
//     //Force GUI options to reflect design pathways
//     if (hasInput) {
//       portDesignOutput.show();
//       //We can't do WtoW design if input is included
//       portWtoWDesignOutput.hide();
//       if (portDesignOpt.getValue() == 4 || portDesignOpt.getValue() == 5) {
//         portDesignOpt.setValue(0,0);
//         theMsg->printf("%s : Cannot do a Manifold-to-Manifold design with input!",getName());
//       }
//       //Change the tdScene with option change
//       if (portDesignOpt.isNew()) renderTopoDesignScene();
//     } else {
//       //Set up WtoW design options
//       portDesignOutput.hide();
//       portWtoWDesignOutput.show();
//     }
//
//
//
//     //------------------------------------------------------------------------------------------------------
//     //Functions for manipulating data:
//     //------------------------------------------------------------------------------------------------------
//     int dataOrbitID = funcOrbitID.getValue();
//     int dataManID = funcManifoldID.getValue();
//     //CLEAR the Manifolds of an orbit
//     if(clearFuncs.wasHit(0)) {
//         int numManifolds = mData->getNumManifolds();
//         bool newData = false;
//         for(int mID=0;mID<numManifolds;mID++) {
//             int baseOrbitID = mData->theManifoldData.mapManifolds[mID].fpdOrbitIdx;
//             if(baseOrbitID == dataOrbitID) {
//                 //If data is there, scrub it...
//                 if(mData->theManifoldData.mapManifolds[mID].isComplete) {
//                     mData->theManifoldData.mapManifolds[mID].clearManifold();
//                     newData = true;
//                 }
//             }
//         }
//         //Render is called since data is new.
//         if(newData) mData->touch();
//     }
//
//     //CLEAR the Manifold of a given ID
//     if(clearFuncs.wasHit(1)) {
//         //Check if manifold actually exists
//         int numManifolds = mData->getNumManifolds();
//         if(dataManID>=0 && dataManID<numManifolds) {
//             //If it does exist, clear that specific manifold
//             mData->theManifoldData.mapManifolds[dataManID].clearManifold();
//             //Indicate that the data is modified in Avizo:
//             mData->touch();
//         } else {
//             //Prompt user and do nothing
//             theMsg->printf("%s: Unable to Clear Manifold ID = %d. Does Not Exist!",__FILE__,dataManID);
//         }
//
//         //Render is called since data is new.
//     }
//
//     //EXTRACT the manifolds of the selected orbit index int new ManifoldData
//     if(extractFuncs.wasHit(0)) {
//         //Create object (blank, missing fpdata)
//         ManifoldData *newMData = new ManifoldData;
//         //Assign the fixed point data object (just the same)
//         newMData->portFPDataFile.setFilename(mData->portFPDataFile.getFilename());
//         newMData->theManifoldData.setFixedPointData(fpData->theData);
//         //Pull out just the Single Orbit 's manifolds
//         int numManifolds = mData->getNumManifolds();
//         int newManID = 0;
//         for(int mID=0;mID<numManifolds;mID++) {
//             int baseOrbitID = mData->theManifoldData.mapManifolds[mID].fpdOrbitIdx;
//             if(baseOrbitID == dataOrbitID) {
//               //Create the new data in the std::map<>
//               newMData->theManifoldData.mapManifolds.insert(
//                 std::pair<int,ManifoldData::ManifoldType>
//                 (newManID, ManifoldData::ManifoldType(mData->theManifoldData.mapManifolds[mID]))
//               );
//               //Make some alterations based on subset
//               ManifoldData::ManifoldType& newManifold = newMData->theManifoldData.mapManifolds[newManID];
//               newManifold.manifoldID = newManID;
//               //We keep the same FPData object (so indexing doesn't change)
//               newManID++;
//             }
//         }
//
//         //Register
//         HxData::registerData(newMData, "ExtractedSingleOrbitMan.im");
//         newMData->fpxDataConnection.connect( fpData );
//         newMData->touch();
//         newMData->compute();
//     }
//     //ADD the given orbit ID to list
//     if(extractFuncs.wasHit(1)) {
//         orbitIDList.push_back( dataOrbitID );
//         printListToPort();
//     }
//     //REMOVE the given orbit ID from list
//     if(extractFuncs.wasHit(2)) {
//         //Find the correct OrbitID and erase
//         std::list<int>::iterator lit;
//         for(lit=orbitIDList.begin();lit!=orbitIDList.end();lit++) {
//           if((*lit) == dataOrbitID) {
//             lit = orbitIDList.erase(lit);
//             --lit;
//           }
//         }
//         printListToPort();
//     }
//     //EXTRACT manifolds of indicated orbit list into new ManifoldData
//     if(extractFuncs.wasHit(3)) {
//         //Create object (blank, missing fpdata)
//         ManifoldData *newMData = new ManifoldData;
//         //Assign the fixed point data object (just the same)
//         newMData->portFPDataFile.setFilename(mData->portFPDataFile.getFilename());
//         newMData->theManifoldData.setFixedPointData(fpData->theData);
//         //Pull out just the Single Orbit 's manifolds
//         int numManifolds = mData->getNumManifolds();
//         int newManID = 0;
//         for(int mID=0;mID<numManifolds;mID++) {
//             int baseOrbitID = mData->theManifoldData.mapManifolds[mID].fpdOrbitIdx;
//             //Check the whole orbit list
//             std::list<int>::iterator lit = orbitIDList.begin();
//             bool inList = false;
//             for(;lit!=orbitIDList.end();lit++) {
//               if ( baseOrbitID == (*lit) ) inList=true;
//             }
//
//             //Add manifold to extracted data if in the list
//             if(inList) {
//               //Create the new data in the std::map<>
//               newMData->theManifoldData.mapManifolds.insert(
//                 std::pair<int,ManifoldData::ManifoldType>
//                 (newManID, ManifoldData::ManifoldType(mData->theManifoldData.mapManifolds[mID]))
//               );
//               //Make some alterations based on subset
//               ManifoldData::ManifoldType& newManifold = newMData->theManifoldData.mapManifolds[newManID];
//               newManifold.manifoldID = newManID;
//               //We keep the same FPData object (so indexing doesn't change)
//               newManID++;
//             }
//         }
//
//         //Register
//         HxData::registerData(newMData, "ExtractedOrbitListMan.im");
//         newMData->fpxDataConnection.connect( fpData );
//         newMData->touch();
//         newMData->compute();
//
//     }
//
//
//     //-----------------------------------------------------------------------------------------------
//     //Rendering updates when ports are changed
//     //-----------------------------------------------------------------------------------------------
//     if (portData.isNew() ||
//         portXbounds.isNew() || portYbounds.isNew() ||
//         horzMenu.isNew() || vertMenu.isNew() || horzScale.isNew() || vertScale.isNew() ||
//         displayOptions.isNew() || eigenLength.isNew() ||
//         colors.isNew() || stableColormap.isNew() || unstableColormap.isNew() ||
//         portLineWidth.isNew() || portPointSize.isNew() ||
//         filterOptions.isNew() || filterType.isNew() || filterOrbitID.isNew() || filterFixedPoint.isNew() || hcOutput.isItemNew(1)
//        ) {
//       renderManifolds();
//     }
//
//     //-----------------------------------------------------------------------------------------------
//     //Show a selected (and propagated) arc
//     //-----------------------------------------------------------------------------------------------
//     if ((simDirection.isNew() || filterOptions.isNew() || filterType.isNew() ||
//          filterOrbitID.isNew() || filterFixedPoint.isNew() ||
//          returnDisplayOpts.isNew() || returnMarkerScale.isNew() ||
//          simOptions.getValue(0)==1||simOptions.getValue(1)==1||simOptions.getValue(2)==1) ) {
//       //Have to propagate if simulation parameters have changed
//       if (show && simDirection.isNew() ) integrate = true;
//       if (show && simOptions.isNew() && simOptions.getValue(0)==1) integrate = true;
//       if (show && returnDisplayOpts.isNew() && returnDisplayOpts.getValue(2)==1) integrate = true;
//       //If integration triggered, call integration function...
//
//       //Function that propagates and displays trajectory scene
//       arcSimulation();
//     }
//
//     //-----------------------------------------------------------------------------------------------
//     // Rendering changes to TopoDesign scenes :
//     //-----------------------------------------------------------------------------------------------
//     if (hasInput || initStateConnection.isNew() ) {
//       //If a change in rendering options for a valid design, render
//       if (initStateConnection.isNew() ||
//           portYbounds.isNew() || portDesignDisplayOpts.isNew() || portDeltaVUnits.isNew() ||
//           designColors.isNew() || portDesignFont.isNew() ) {
//         renderTopoDesignScene();
//       }
//     }
//     if (!hasInput && augHConn.isValid) {
//       //If a change in rendering options for a valid design, render
//       if (portYbounds.isNew() || portDesignDisplayOpts.isNew() || portDeltaVUnits.isNew() ||
//           designColors.isNew() || portDesignFont.isNew() ) {
//         renderWWTopoDesignScene(false);
//       }
//     }
//
//
//     //-----------------------------------------------------------------------------------------------
//     //Data output buttons for arc information
//     //-----------------------------------------------------------------------------------------------
//     if (simData.wasHit(0)) { //LineSet output
//       //If data exists
//       if(show) {
//         HxLineSet *lineset = new HxLineSet;
//         exportArcToLineSet(lineset);
//         lineset->touchMinMax();
//         lineset->setLabel("ManifoldArc.hx");
//         theObjectPool->addObject(lineset);
//       } else {
//         theMsg->printf("No available data to construct an HxLineSet object.");
//         return;
//       }
//     }
//     if (simData.wasHit(1)) { //Data output (text)
//       theMsg->printf("%s: Data text output under construction!",getName());
//     }
//
//     if (simData.wasHit(2)) { //Output Trajectory object for design
//       //If data exists
//       if(show) {
//         //Run the simulation if not available
//         if(simOptions.getValue(0)!=1) {
//           integrate = true;
//           //Function for propagation
//           simOptions.setValue(0,1);
//           arcSimulation();
//           simOptions.setValue(0,0);
//         }
//         HxTrajectory *arc = new HxTrajectory;
//         exportArcToTrajectory(arc);
//         arc->setLabel("ManifoldArc");
//         theObjectPool->addObject(arc);
//       } else {
//         theMsg->printf("No available data to construct a Trajectory object.");
//         return;
//       }
//     }
//
//     if (simData.wasHit(3)) { //Output Returns as HxCluster
//         //If data exists
//         if(show) {
//             //Run simulation if not available
//             if(simOptions.getValue(0)!=1) {
//                 integrate = true;
//                 simOptions.setValue(0,1);
//                 arcSimulation();
//                 simOptions.setValue(0,0);
//             }
//             HxCluster *rets = new HxCluster;
//             exportReturnsToCluster(rets);
//             rets->setLabel("ManifoldArcReturns");
//             theObjectPool->addObject(rets);
//         } else {
//             theMsg->printf("No available data to construct returns.");
//             return;
//         }
//     }
//
//     //Output Topology Design Maneuver from Input
//     if( portDesignOutput.wasHit(1) ) {
//         //Need to have input
//         if(!hasInput) {
//           theMsg->printf("%s : No InitialState input detected, no maneuver to process.",__FILE__);
//           return;
//         }
//         //Need to have a selection
//         if(!topoDesignSelected) {
//           theMsg->printf("%s : No manifold selection to form topology-based deltaV.",__FILE__);
//           return;
//         }
//
//         //Use input
//         HxInitialState* iState = (HxInitialState*) portData.source();
//         State &theState = iState->data.ic;
//         //Also need the manifold data we picked to formulate name
//         LineSegPair pickedManifold( selectedLineID, selectedPartID );
//         //LineSegPair hcPartnerManifold( hcPartnerLineID, hcPartnerPartID );
//         //LineID,SegID -> ManifoldID
//         IntPair manSegPair;
//         //IntPair hcManSegPair;
//         if(selectedUnstable) {
//           //Get (manifoldID,segmentID)
//           manSegPair = uLinesetToManifold[pickedManifold];
//         } else {
//           manSegPair = sLinesetToManifold[pickedManifold];
//         }
//         /*//Heteroclinic piece - Ignore
//         if(hcOutput.getValue(0)==1) {
//           if(selectedUnstable) {
//             //Get (manifoldID,segmentID) for HC partner
//             hcManSegPair = sLinesetToManifold[hcPartnerManifold];
//           } else {
//             hcManSegPair = uLinesetToManifold[hcPartnerManifold];
//           }
//         }*/
//
//         //Get the Manifold and Segment
//         int mID = manSegPair.first; //int segID = manSegPair.second;
//         const ManifoldData::ManifoldType& theManifold = mData->theManifoldData.mapManifolds[mID];
//         //const ManifoldData::ManifoldSeg& theSeg = theManifold.segments[segID];
//         int fpOrbitID = theManifold.fpdOrbitIdx;
//         //Fill out Maneuver object
//         HxManeuver *dvObj = new HxManeuver;
//         dvObj->data.dv = theState; //For location
//         dvObj->data.dv.xd = deltaV[0];
//         dvObj->data.dv.yd = deltaV[1];
//         dvObj->data.dv.zd = deltaV[2];
//         QString oIDText,nameStr;
//         oIDText.setNum(fpOrbitID);
//         switch(portDesignOpt.getValue()) {
//           case 0 :
//             dvObj->setLabel("MapDV_StateToWs.dv");
//             nameStr = "Entrance onto $W^s$ of Orbit " + oIDText + " from State";
//             break;
//           case 1 :
//             dvObj->setLabel("MapDV_WuToState.dv");
//             nameStr = "Departure from $W^u$ of Orbit " + oIDText + " to State";
//             break;
//           case 2 :
//             dvObj->setLabel("MapDV_StateToWu.dv");
//             //nameStr = "State to $W^u$ of Orbit " + oIDText;
//             nameStr = "Depature from $W^u$ of Orbit " + oIDText + " to State"; //Only (F) design names
//             break;
//           case 3 :
//             dvObj->setLabel("MapDV_WsToState.dv");
//             //nameStr = "$W^s$ of Orbit " + oIDText " to State";
//             nameStr = "Entrance onto $W^s$ of Orbit " + oIDText + " from State"; //Only (F) design names
//             break;
//         }
//         dvObj->data.name = nameStr.toLocal8Bit().data();
//         dvObj->updatePorts();
//         dvObj->touch();
//         theObjectPool->addObject(dvObj);
//         dvObj->compute();
//     }
//
//     //Output an trajectory for the topology-based design from InitialState Input
//     if( portDesignOutput.wasHit(0) ) {
//         //Need to have input
//         if(!hasInput) {
//           theMsg->printf("%s : No InitialState input detected, no topology-based arc to process.",__FILE__);
//           return;
//         }
//         //Need to have a selection
//         if(!topoDesignSelected) {
//           theMsg->printf("%s : No manifold selection to form topology-based deltaV and new arc.",__FILE__);
//           return;
//         }
//
//         //Construct the Trajectory (upstream object only!)
//         if(simOptions.getValue(0)!=1) {
//           integrate = true;
//           //Function for propagation of selection
//           simOptions.setValue(0,1);
//           arcSimulation();
//           simOptions.setValue(0,0);
//         }
//         //Gather the orbit index for name
//         LineSegPair pickedManifold( selectedLineID, selectedPartID );
//         //LineID,SegID -> ManifoldID
//         IntPair manSegPair;
//         //IntPair hcManSegPair;
//         if(selectedUnstable) {
//           //Get (manifoldID,segmentID)
//           manSegPair = uLinesetToManifold[pickedManifold];
//         } else {
//           manSegPair = sLinesetToManifold[pickedManifold];
//         }
//         //Get the Manifold and Segment
//         int mID = manSegPair.first; //int segID = manSegPair.second;
//         const ManifoldData::ManifoldType& theManifold = mData->theManifoldData.mapManifolds[mID];
//         //const ManifoldData::ManifoldSeg& theSeg = theManifold.segments[segID];
//         int fpOrbitID = theManifold.fpdOrbitIdx;
//         int fpPointID = theManifold.fpdPointIdx;
//         HxTrajectory *wArc = new HxTrajectory;
//         exportArcToTrajectory(wArc);
//         //Create name based on design option
//         QString oIDText,fpIDText,nameStr;
//         oIDText.setNum(fpOrbitID);
//         fpIDText.setNum(fpPointID);
//         switch(portDesignOpt.getValue()) {
//           case 0 :
//             wArc->setLabel("WsArc.traj");
//             nameStr = "$W^s$ of Orbit " + oIDText + " (" + fpIDText + ")";
//             break;
//           case 1 :
//             wArc->setLabel("WuArc.traj");
//             nameStr = "$W^u$ of Orbit " + oIDText + " (" + fpIDText + ")";
//             break;
//           case 2 :
//             wArc->setLabel("WuArc.traj");
//             nameStr = "$W^u$ of Orbit " + oIDText + " (" + fpIDText + ")";
//             break;
//           case 3 :
//             wArc->setLabel("WsArc.traj");
//             nameStr = "$W^s$ of Orbit " + oIDText + " (" + fpIDText + ")";
//             break;
//         }
//         wArc->data.name = nameStr.toLocal8Bit().data();
//         wArc->updatePorts();
//         theObjectPool->addObject(wArc);
//     }
//
//     //Output Manifold to Manifold Topology Design Objects
//     if( portWtoWDesignOutput.wasHit(0) ) {
//         //Make sure design was created
//         if(!augHConn.isValid) {
//             theMsg->printf("%s : Need a valid augmented connection computed (manifold to manifold).",__FILE__);
//             return;
//         }
//         //Output objects:
//         // 1) Starting Manifold Arc from periodic orbit A
//         // 2) InitialState for end of initial manifold
//         // 3) Maneuver for map-based translation
//         // 4) InitialState for start of next manifold
//         // 5) Next Manifold Arc back to periodic orbit B
//
//         //Build Names
//         std::vector<QString> nameStrs(5);
//         QString oIDText0, oIDText1;
//         oIDText0.setNum( augHConn.orbitID0 );
//         oIDText1.setNum( augHConn.orbitID1 );
//         QString fpIDText0, fpIDText1;
//         fpIDText0.setNum( augHConn.fpID0 );
//         fpIDText1.setNum( augHConn.fpID1 );
//         QString traj0Label, traj1Label;
//         QString s0Label, s1Label;
//         if(portDesignOpt.getValue() == 4) {
//           traj0Label = "AugHc_WuArc.traj";
//           traj1Label = "AugHc_WsArc.traj";
//           s0Label = "AugHc_WuState.ist";
//           s1Label = "AugHc_WsState.ist";
//           QString ahcStr = "Aug. $\\mathcal{H}_c [" + oIDText0 + "\\mapsto " + oIDText1 + "]$: ";
//           nameStrs[0] = ahcStr + "$W^u$ of Orbit " + oIDText0 + " (" + fpIDText0 + ")";
//           nameStrs[1] = ahcStr + "State at end of $W^u$ of Orbit " + oIDText0 + " (" + fpIDText0 + ")";
//           nameStrs[2] = ahcStr + "Maneuver $W^u$ to $W^s$";
//           nameStrs[3] = ahcStr + "State at start of $W^s$ of Orbit " + oIDText1 + " (" + fpIDText1 + ")";
//           nameStrs[4] = ahcStr + "$W^s$ of Orbit " + oIDText1 + " (" + fpIDText1 + ")";
//         } else {
//           traj1Label = "AugHc_WuArc.traj";
//           traj0Label = "AugHc_WsArc.traj";
//           s1Label = "AugHc_WuState.ist";
//           s0Label = "AugHc_WsState.ist";
//           QString ahcStr = "Aug. $\\mathcal{H}_c [" + oIDText1 + "\\mapsto " + oIDText0 + "]$: ";
//           nameStrs[0] = ahcStr + "$W^s$ of Orbit " + oIDText0 + " (" + fpIDText0 + ")";
//           nameStrs[1] = ahcStr + "State at start of $W^s$ of Orbit " + oIDText0 + " (" + fpIDText0 + ")";
//           nameStrs[2] = ahcStr + "Maneuver $W^u$ to $W^s$";
//           nameStrs[3] = ahcStr + "State at end of $W^u$ of Orbit " + oIDText1 + " (" + fpIDText1 + ")";
//           nameStrs[4] = ahcStr + "$W^u$ of Orbit " + oIDText1 + " (" + fpIDText1 + ")";
//         }
//
//         //Check the simulation option to force propagations
//         bool simOptOff = false;
//         if(simOptions.getValue(0)!=1) {
//           simOptOff = true;
//           simOptions.setValue(0,1);
//         }
//         //Need to integrate to construct the trajectories
//         integrate = true;
//
//         //-------------------------------
//         // 1) Starting Manifold arc:
//         //-------------------------------
//         //Modify member to trigger integration of correct arcs
//         selectedUnstable = augHConn.fwd;
//         selectedLineID = augHConn.man0LineIdx;
//         selectedPartID = augHConn.man0PartIdx;
//         selectedLinearParameter = augHConn.tau0;
//         //Simulate - with Option to skip computation of linearParam and use set value
//         arcSimulation(true);
//         double connectionTOF = 0.0;
//         HxTrajectory *arc0 = new HxTrajectory;
//         exportArcToTrajectory(arc0);
//         arc0->setLabel(traj0Label.toLocal8Bit().data());
//         arc0->data.name = nameStrs[0].toLocal8Bit().data();
//         arc0->updatePorts();
//         theObjectPool->addObject( arc0 );
//         connectionTOF += std::fabs( arc0->data.getTimeOfFlight() );
//
//         //----------------------------------------
//         // 2) State object at end of starting arc
//         //----------------------------------------
//         HxInitialState *s0 = new HxInitialState;
//         s0->data.name =  nameStrs[1].toLocal8Bit().data();
//         s0->setData( augHConn.manifoldState0 );
//         s0->setLabel( s0Label.toLocal8Bit().data() );
//         theObjectPool->addObject( s0 );
//
//         //----------------------------------------
//         // 3) Maneuver
//         //----------------------------------------
//         HxManeuver *dvObj = new HxManeuver;
//         dvObj->data.dv = augHConn.deltaV;
//         dvObj->data.name =  nameStrs[2].toLocal8Bit().data();
//         dvObj->updatePorts();
//         dvObj->setLabel( "AugHc_DeltaV.dv" );
//         theObjectPool->addObject( dvObj );
//
//         //----------------------------------------
//         // 4) State object at start of last arc
//         //----------------------------------------
//         HxInitialState *s1 = new HxInitialState;
//         s1->data.name =  nameStrs[3].toLocal8Bit().data();
//         s1->setData( augHConn.manifoldState1 );
//         s1->setLabel( s1Label.toLocal8Bit().data() );
//         theObjectPool->addObject( s1 );
//
//         //-------------------------------
//         // 5) Last Manifold arc:
//         //-------------------------------
//         //Modify member to trigger integration of correct arcs
//         selectedUnstable = (augHConn.fwd)? false : true; //Propagate from the other set
//         selectedLineID = augHConn.man1LineIdx;
//         selectedPartID = augHConn.man1PartIdx;
//         selectedLinearParameter = augHConn.tau1;
//         //Simulate - with Option to skip computation of linearParam and use set value
//         integrate = true;
//         arcSimulation(true);
//         HxTrajectory *arc1 = new HxTrajectory;
//         exportArcToTrajectory(arc1);
//         arc1->setLabel(traj1Label.toLocal8Bit().data());
//         arc1->data.name = nameStrs[4].toLocal8Bit().data();
//         arc1->updatePorts();
//         theObjectPool->addObject( arc1 );
//         connectionTOF += std::fabs( arc1->data.getTimeOfFlight() );
//
//         //Modify the resulting tof as it is a better estimate
//         augHConn.tof = connectionTOF;
//         //Done with this
//         integrate = false;
//         //Prompt user
//         theMsg->printf("%s: Augmented Connection:  ToF = %.15f (nd)  DeltaV = %.15f (nd)",getName(),augHConn.tof,augHConn.deltaV.vmag());
//         //Reset switch if necessary
//         if(simOptOff) {
//           simOptions.setValue(0,0);
//           //Call again to remove
//           arcSimulation();
//         }
//     }
//     if( portWtoWDesignOutput.wasHit(1) ) {
//         //Make sure design was created
//         if(!augHConn.isValid) {
//             theMsg->printf("%s : Need a valid augmented connection computed (manifold to manifold).",__FILE__);
//             return;
//         }
//         //Print WtoWDesign object to screen
//         augHConn.printToConsole();
//     }
//     if( portWtoWDesignOutput.wasHit(2) ) {
//         //Clear the topology design scene
//         renderWWTopoDesignScene(true); //Remove's from display
//         //Clear the WtoWDesign
//         augHConn.reset();
//     }
//
//
//     //Heteroclinic/Homoclinic connection extraction
//     if (hcOutput.getValue(2) == 1) { //Button was hit
//         //Make sure a connection was detected
//         if (hcOutput.getValue(0)) {
//             bool simOptOff = false;
//             if(simOptions.getValue(0)!=1) {
//               simOptOff = true;
//               simOptions.setValue(0,1);
//             }
//             //Simulate unstable
//             integrate = true;
//             //theMsg->printf("HC: Simulating unstable arc...");
//             arcSimulation();
//
//             //Time of flight for connection
//             double connectionTOF = 0.0;
//             //Output
//             HxTrajectory *uArc = new HxTrajectory;
//             exportArcToTrajectory(uArc);
//             uArc->setLabel("HC_UnstableArc");
//             //Construct name
//             QString oIDText,fpIDText,nameStr;
//             oIDText.setNum(hhConn.orbitID0);
//             fpIDText.setNum(hhConn.fpID0);
//             nameStr = "$W^u$ of Orbit " + oIDText + " (" + fpIDText + ") for $\\mathcal{H}_c$";
//             uArc->data.name = nameStr.toLocal8Bit().data();
//             uArc->updatePorts();
//             theObjectPool->addObject(uArc);
//             connectionTOF += std::fabs( uArc->data.getTimeOfFlight() );
//             //Switch to stable & simulate
//             selectedUnstable = false;
//             int tempLine = selectedLineID;
//             int tempPart = selectedPartID;
//             selectedLineID = hcPartnerLineID;
//             selectedPartID = hcPartnerPartID;
//             hcPartnerLineID = tempLine;
//             hcPartnerPartID = tempPart;
//             integrate = true;
//             //theMsg->printf("HC: Simulating stable arc...");
//             arcSimulation();
//             //Output
//             HxTrajectory *sArc = new HxTrajectory;
//             exportArcToTrajectory(sArc);
//             sArc->setLabel("HC_StableArc");
//             oIDText.setNum(hhConn.orbitID1);
//             fpIDText.setNum(hhConn.fpID1);
//             nameStr = "$W^s$ of Orbit " + oIDText + " (" + fpIDText + ") for $\\mathcal{H}_c$";
//             sArc->data.name = nameStr.toLocal8Bit().data();
//             sArc->updatePorts();
//             theObjectPool->addObject(sArc);
//             connectionTOF += std::fabs( sArc->data.getTimeOfFlight() );
//             //This is now more accurate information, so store in object
//             hhConn.tof = connectionTOF;
//             //Done with this
//             integrate = false;
//             //Prompt the user about the transfer ToF
//             theMsg->printf("%s: Heteroclinic connection time-of-flight = %.15f (non-dim)",getName(),connectionTOF);
//             //Reset switch if necessary
//             if(simOptOff) {
//               simOptions.setValue(0,0);
//               //Call again (likely to remove arc)
//               arcSimulation();
//             }
//         } else {
//             theMsg->printf("%s: No Heteroclinic/Homoclinic Connection detected",getName());
//             return;
//         }
//     }
//     //Heteroclinic/Homoclinic connection data output
//     if (hcOutput.getValue(3) == 1) { //Button was hit
//         //Make sure it's a valid connection point
//         if (hcOutput.getValue(0)) {
//           //Print WWDesign for hhConn
//           hhConn.printToConsole();
//         } else {
//           theMsg->printf("%s: No Heteroclinic/Homoclinic Connection detected",getName());
//           return;
//         }
//     }
//
// }
//
// /// Perform the propagation step, render, and show
// void DisplayManifolds::arcSimulation(const bool useLinearParam)
// {
//     //Check for input data
//     ManifoldData *mData = (ManifoldData*) portData.source();
//     if(!mData) {
//         theMsg->printf("%s: Cannot render without valid manifold data",__FILE__);
//         return;
//     }
//
//     std::vector<State> uItStates, dItStates;
//     if (integrate) {
//         //Run the integration
//         uStates.clear(); uTimes.clear();
//         dStates.clear(); dTimes.clear();
//         //Call the propagator (stores both full states and iterates)
//         integrateManifoldArc(useLinearParam);
//         integrate = false;
//     }
//
//     //Gather manifold information if necessary
//     LineSegPair pickedManifold( selectedLineID, selectedPartID );
//     //LineID,SegID -> ManifoldID
//     IntPair manSegPair;
//     if(selectedUnstable) {
//         //Get (manifoldID,segmentID)
//         manSegPair = uLinesetToManifold[pickedManifold];
//     } else {
//         manSegPair = sLinesetToManifold[pickedManifold];
//     }
//     int mID = manSegPair.first; //int segID = manSegPair.second;
//
//     //Clear the current scene
//     trajScene->removeAllChildren();
//
//
//     //Render the trajectory
//     if(show && simOptions.getValue(0)==1) {
//         //Render with SoLineSet
//         SoSeparator *lineSep = new SoSeparator;
//         SoLineSet *line = new SoLineSet;
//         SoLineSet *dline = new SoLineSet;
//         SoCoordinate3 *lineCoords = new SoCoordinate3;
//         SoCoordinate3 *dlineCoords = new SoCoordinate3;
//         SoDrawStyle *lStyle = new SoDrawStyle;
//         SoDrawStyle *dlStyle = new SoDrawStyle;
//         lStyle->lineWidth = portLineWidth.getValue();
//         dlStyle->lineWidth = portLineWidth.getValue();
//         //Material
//         SoMaterial *lMat = new SoMaterial;
//         SoMaterial *dlMat = new SoMaterial;
//         McColor *lCol = new McColor;
//         //Color selected by which manifold we are on (Wu or Ws)
//         int fOrbitID = mData->theManifoldData.mapManifolds[mID].fpdOrbitIdx;
//         if(selectedUnstable) {
//           //colors.getColor(*lCol,2);
//           //Colormap
//           SbColor c = unstableColormap.getColor(fOrbitID);
//           lCol->setValue(c[0],c[1],c[2]);
//         } else {
//           //colors.getColor(*lCol,1);
//           //Colormap
//           SbColor c = stableColormap.getColor(fOrbitID);
//           lCol->setValue(c[0],c[1],c[2]);
//         }
//
//         //Compute compliment
//         float h,s,v; h = s = v = 0.0;
//         lCol->getHSVValue(h,s,v);
//         h += 0.5; //Add half hue (keep s,v the same)
//         if (h>1.0) h-=1.0;
//         McColor *complCol = new McColor;
//         complCol->setHSVValue(h,s,v);
//
//         //Assign color of trajectory objects based on selection
//         switch (simDirection.getValue() ) {
//             case 1 :
//                 //Upstream only
//                 lMat->diffuseColor.setValue((*lCol)[0],(*lCol)[1],(*lCol)[2]);
//                 if (simRevColor.getValue() == 1) {
//                     dlMat->diffuseColor.setValue((*complCol)[0],(*complCol)[1],(*complCol)[2]);
//                 } else {
//                     dlMat->diffuseColor.setValue((*lCol)[0],(*lCol)[1],(*lCol)[2]);
//                 }
//                 break;
//             default :
//                 //Both - Make Downstream the original color from color list
//                 dlMat->diffuseColor.setValue((*lCol)[0],(*lCol)[1],(*lCol)[2]);
//                 if (simRevColor.getValue() == 1) {
//                     lMat->diffuseColor.setValue((*complCol)[0],(*complCol)[1],(*complCol)[2]);
//                 } else {
//                     lMat->diffuseColor.setValue((*lCol)[0],(*lCol)[1],(*lCol)[2]);
//                 }
//                 break;
//         }
//
//         //Add objects to manifold display
//         lineSep->addChild(lStyle);
//         lineSep->addChild(lMat);
//         lineSep->addChild(lineCoords);
//         lineSep->addChild(line);
//
//         //Set the UPSTREAM Line coords
//         int numStates = (int) uStates.size();
//         for(int ii=0;ii<numStates;ii++) {
//           lineCoords->point.set1Value(ii,uStates[ii].x,uStates[ii].y,uStates[ii].z);
//         }
//         line->numVertices.set1Value(0,numStates);
//
//         //Set the DOWNSTREAM data
//         lineSep->addChild(dlStyle);
//         lineSep->addChild(dlMat);
//         lineSep->addChild(dlineCoords);
//         lineSep->addChild(dline);
//         //Set the DOWNSTREAM Line coords
//         numStates = (int) dStates.size();
//         for(int ii=0;ii<numStates;ii++) {
//           dlineCoords->point.set1Value(ii,dStates[ii].x,dStates[ii].y,dStates[ii].z);
//         }
//         dline->numVertices.set1Value(0,numStates);
//
//         //Show
//         trajScene->addChild(lineSep);
//
//     }
//
//     //Render the Progenitor State
//     /*if(show && simOptions.getValue(3)==1) {
//         McColor *wCol = new McColor;
//         //Color selected by which manifold we are on (Wu or Ws)
//         if(selectedUnstable) {
//           colors.getColor(*wCol,2);
//         } else {
//           colors.getColor(*wCol,1);
//         }
//         //Compute compliment
//         float h,s,v; h = s = v = 0.0;
//         wCol->getHSVValue(h,s,v);
//         h += 0.5; //Add half hue (keep s,v the same)
//         if (h>1.0) h-=1.0;
//         McColor *compwCol = new McColor;
//         compwCol->setHSVValue(h,s,v);
//
//         //Show the Progenitor State
//         SoSeparator *psSep = new SoSeparator;
//         SbVec3f *psMarkerCoord = new SbVec3f(progStatePoint);
//         SbColor *psColor = new SbColor( (*wCol)[0],(*wCol)[1],(*wCol)[2] );
//         //Switch color if UPSTREAM is compliment
//         if(simRevColor.getValue() == 1) psColor->setHSVValue(h,s,v);
//         int psMarkerIdx[1] = { SoMarkerSet::DIAMOND_FILLED_5_5 };
//         SoCoordinate3 *psCoord = new SoCoordinate3;
//         psCoord->point.setValues(0,1,psMarkerCoord);
//         SoMarkerSet *psMarkerSet = new SoMarkerSet;
//         psMarkerSet->markerIndex.setValues(0,1,psMarkerIdx);
//         psMarkerSet->markerScale.setValue( portMarkerScale.getValue() );
//         SoMaterialBinding *psMatBind = new SoMaterialBinding;
//         psMatBind->value = SoMaterialBinding::PER_VERTEX;
//         SoMaterial *psMarkerMat = new SoMaterial;
//         psMarkerMat->diffuseColor.setValues(0,1,psColor);
//         psSep->addChild(psCoord);
//         psSep->addChild(psMatBind);
//         psSep->addChild(psMarkerMat);
//         psSep->addChild(psMarkerSet);
//         //Add Progenitor State to trajectory scene
//         trajScene->addChild(psSep);
//     }*/
//
//
//     //Render the iterations
//     if(show && simOptions.getValue(1)==1) {
//
//         //Create rendering objects
//         SoSeparator *itSep = new SoSeparator;
//         SoDepthOffset *itDepthOffset = new SoDepthOffset();
//         itSep->addChild(itDepthOffset);
//
//         McColor *wCol = new McColor;
//         //Color selected by which manifold we are on (Wu or Ws)
//         int fOrbitID = mData->theManifoldData.mapManifolds[mID].fpdOrbitIdx;
//         if(selectedUnstable) {
//           //Color button
//           //colors.getColor(*wCol,2);
//           //Colormap
//           SbColor c = unstableColormap.getColor(fOrbitID);
//           wCol->setValue(c[0],c[1],c[2]);
//         } else {
//           //Color button
//           //colors.getColor(*wCol,1);
//           //Colormap
//           SbColor c = stableColormap.getColor(fOrbitID);
//           wCol->setValue(c[0],c[1],c[2]);
//         }
//         //Compute compliment
//         float h,s,v; h = s = v = 0.0;
//         wCol->getHSVValue(h,s,v);
//         h += 0.5; //Add half hue (keep s,v the same)
//         if (h>1.0) h-=1.0;
//         McColor *compwCol = new McColor;
//         compwCol->setHSVValue(h,s,v);
//
//         //Show the Selection Point
//         //---------------------------------------------------------------------------------------
//         if (returnDisplayOpts.getValue(0)==1) {
//             SbVec3f *ppMarkerCoord = new SbVec3f(pickedPoint);
//             SbColor *markerColor = new SbColor((*wCol)[0],(*wCol)[1],(*wCol)[2]);
//             //Change color if connection is detected
//             if(hcOutput.getValue(0)==1) {
//                 McColor hcCol;
//                 hcOutput.getColor(1,hcCol);
//                 markerColor->setValue(hcCol[0],hcCol[1],hcCol[2]);
//             }
//             int markerIndex[1] = { SoMarkerSet::SATELLITE_FILLED_9_9};
//             SoCoordinate3 *markerCoord = new SoCoordinate3;
//             markerCoord->point.setValues(0,1,ppMarkerCoord);
//             SoMarkerSet *markerSet = new SoMarkerSet;
//             markerSet->markerIndex.setValues(0,1,markerIndex);
//             markerSet->markerScale.setValue( returnMarkerScale.getValue() );
//             SoMaterialBinding *matBind = new SoMaterialBinding;
//             matBind->value = SoMaterialBinding::PER_VERTEX;
//             SoMaterial *markerMat = new SoMaterial;
//             markerMat->diffuseColor.setValues(0,1,markerColor);
//             itSep->addChild(markerCoord);
//             itSep->addChild(matBind);
//             itSep->addChild(markerMat);
//             itSep->addChild(markerSet);
//
//             //For Hc's, render the time of flight text (estimate)
//             if (hcOutput.getValue(0)==1) {
//                 SoMaterial *tofMat = new SoMaterial;
//                 if(hhConn.isValid) {
//                     McColor fontCol( portDesignFont.getFontColor() );
//                     tofMat->diffuseColor.setValue(fontCol.r,fontCol.g,fontCol.b);
//                 } else {
//                     //Indicate an invalid point by red!
//                     tofMat->diffuseColor.setValue(SbColor(1.0,0.0,0.0));
//                 }
//                 //Font
//                 SoFont *textFont = new SoFont;
//                 QString familyName = portDesignFont.getFontName();
//                 bool isBold = portDesignFont.isBoldFont();
//                 bool isItalic = portDesignFont.isItalicFont();
//                 if (isBold && isItalic) {
//                   familyName.append(" : Bold Italic");
//                 } else if (isBold) {
//                   familyName.append(" : Bold");
//                 } else if (isItalic) {
//                   familyName.append(" : Italic");
//                 }
//                 QByteArray fbarray = familyName.toLocal8Bit();
//                 textFont->name.setValue(fbarray.data());
//                 textFont->size.setValue(portDesignFont.getFontSize());
//                 //Add to sep to control all subsequent SoNodes
//                 SoSeparator *textSep = new SoSeparator;
//                 //textSep->addChild(textDOffset);
//                 textSep->addChild(tofMat);
//                 textSep->addChild(textFont);
//                 //Indicate a translation to define the location
//                 SoTransform *textXform = new SoTransform;
//                 textXform->translation.setValue( ppMarkerCoord[0] - SbVec3f(0,0,0) );
//                 //The actual Text
//                 QString tofTextStr;
//                 double value = 0.0;
//                 if (hhConn.isValid) {
//                     //Convert the approx time of flight in days
//                     value = sys->tstar * hhConn.tof; //days
//                     tofTextStr.setNum(value,'g',5);
//                     tofTextStr += " days";
//
//                 } else {
//                     tofTextStr = "N/A";
//                 }
//                 //Move off the point a bit and add translation with spaces:
//                 QString fullText;
//                 //if (portDesignDisplayOpts.getValue(1)==1) {
//                   fullText = "  " + QString(QChar(0x0394)) + "T= " + tofTextStr;
//                 //}
//                 SoText2 *tofText = new SoText2;
//                 SbString soStr; soStr.fromUtf16(fullText.utf16());
//                 tofText->string.setValue(soStr);
//                 //Add to scene
//                 textSep->addChild(textXform);
//                 textSep->addChild(tofText);
//                 itSep->addChild(textSep);
//             } //End Hc time of flight display
//         }
//
//         //Build integration objects for reference
//         //----------------------------------------------------------------------------------------
//         rhs_type theEOMs(C,mup);  //The right-hand side:  CR3BP EoMs
//         section_type theSection(theEOMs); //Hyperplane: y=0 by planar_section
//
//
//         //Store objects as map-space coords, so we need the map-space to world transform
//         //Map Projection:  Defines a Transform Node that applies to all rendered objects
//         //---------------------------------------------------------------------------------------
//         SoTransform* theXForm = new SoTransform();
//         mapSpaceProjector.compute(horzMenu,vertMenu,horzScale,vertScale,theXForm);
//         itSep->addChild(theXForm);
//
//
//         //Show the Iterates (and sub-iterates) if selected
//         //---------------------------------------------------------------------------------------
//         if (returnDisplayOpts.getValue(1)==1) {
//
//           //UPSTREAM Iterates first (Marked with Circles)
//           //-------------------------------------------------------------------------------------
//           int NUM_ITERS = (int) uIterates.size();
//           SbVec3f *iterCoords = new SbVec3f[NUM_ITERS];
//           SbColor *itersColor = new SbColor[NUM_ITERS];
//           int *iterMarkerIndex = new int[NUM_ITERS];
//           for (int i=0; i<NUM_ITERS; i++) {
//               nvis::vec2 x = uIterates[i];
//               SbVec3f mapPoint(x[0],x[1],0.0);
//               iterCoords[i].setValue(mapPoint[0],mapPoint[1],mapPoint[2]);
//               //Setting the color  - DOWNSTREAM is always the main color
//               if (simRevColor.getValue() == 1 ) {
//                 itersColor[i].setValue((*compwCol)[0],(*compwCol)[1],(*compwCol)[2]);
//               } else {
//                 itersColor[i].setValue((*wCol)[0],(*wCol)[1],(*wCol)[2]);
//               }
//               iterMarkerIndex[i] = SoMarkerSet::CIRCLE_FILLED_5_5;
//           }
//           SoCoordinate3 *iterCoord3 = new SoCoordinate3;
//           iterCoord3->point.setValues(0,NUM_ITERS,iterCoords);
//           SoMarkerSet *iterMarkerSet = new SoMarkerSet;
//           iterMarkerSet->markerIndex.setValues(0,NUM_ITERS,iterMarkerIndex);
//           iterMarkerSet->markerScale.setValue( returnMarkerScale.getValue() );
//           SoMaterialBinding *iterMatBind = new SoMaterialBinding;
//           iterMatBind->value = SoMaterialBinding::PER_VERTEX;
//           SoMaterial *iterMarkerMat = new SoMaterial;
//           iterMarkerMat->diffuseColor.setValues(0,NUM_ITERS,itersColor);
//           //Add UPSTREAM to scene
//           itSep->addChild(iterCoord3);
//           itSep->addChild(iterMatBind);
//           itSep->addChild(iterMarkerMat);
//           itSep->addChild(iterMarkerSet);
//
//           //DOWNSTREAM Iterates (Marked with Diamonds)
//           //-------------------------------------------------------------------------------------
//           NUM_ITERS = (int) dIterates.size();
//           SbVec3f *diterCoords = new SbVec3f[NUM_ITERS];
//           SbColor *ditersColor = new SbColor[NUM_ITERS];
//           int *diterMarkerIndex = new int[NUM_ITERS];
//           for (int i=0; i<NUM_ITERS; i++) {
//               nvis::vec2 x = dIterates[i];
//               SbVec3f mapPoint(x[0],x[1],0.0);
//               diterCoords[i].setValue(mapPoint[0],mapPoint[1],mapPoint[2]);
//               //Setting the color  - DOWNSTREAM is always the main color
//               ditersColor[i].setValue((*wCol)[0],(*wCol)[1],(*wCol)[2]);
//               diterMarkerIndex[i] = SoMarkerSet::DIAMOND_FILLED_5_5;
//           }
//           SoCoordinate3 *diterCoord3 = new SoCoordinate3;
//           diterCoord3->point.setValues(0,NUM_ITERS,diterCoords);
//           SoMarkerSet *diterMarkerSet = new SoMarkerSet;
//           diterMarkerSet->markerIndex.setValues(0,NUM_ITERS,diterMarkerIndex);
//           diterMarkerSet->markerScale.setValue( returnMarkerScale.getValue() );
//           SoMaterialBinding *diterMatBind = new SoMaterialBinding;
//           diterMatBind->value = SoMaterialBinding::PER_VERTEX;
//           SoMaterial *diterMarkerMat = new SoMaterial;
//           diterMarkerMat->diffuseColor.setValues(0,NUM_ITERS,ditersColor);
//           //Add UPSTREAM to scene
//           itSep->addChild(diterCoord3);
//           itSep->addChild(diterMatBind);
//           itSep->addChild(diterMarkerMat);
//           itSep->addChild(diterMarkerSet);
//
//         }
//
//         //Render the sub-iterates
//         //---------------------------------------------------------------------------------------
//         if (returnDisplayOpts.getValue(2)==1) {
//
//           //UPSTREAM Sub-Iterates first (Marked with Open Circles)
//           //-------------------------------------------------------------------------------------
//           int NUM_ITERS = (int) uSubIts.size();
//           SbVec3f *iterCoords = new SbVec3f[NUM_ITERS];
//           SbColor *itersColor = new SbColor[NUM_ITERS];
//           int *iterMarkerIndex = new int[NUM_ITERS];
//           for (int i=0; i<NUM_ITERS; i++) {
//               vec42 y(0);
//               for(int j=0;j<6;j++) y[j] = uSubIts[i][j];
//               //This may need adjust for different maps
//               nvis::vec2 x = theSection.project(y).first;
//               SbVec3f mapPoint(x[0],x[1],0.0);
//               iterCoords[i].setValue(mapPoint[0],mapPoint[1],mapPoint[2]);
//               //Setting the color  - DOWNSTREAM is always the main color
//               if (simRevColor.getValue() == 1 ) {
//                 itersColor[i].setValue((*compwCol)[0],(*compwCol)[1],(*compwCol)[2]);
//               } else {
//                 itersColor[i].setValue((*wCol)[0],(*wCol)[1],(*wCol)[2]);
//               }
//               iterMarkerIndex[i] = SoMarkerSet::CIRCLE_LINE_5_5;
//           }
//           SoCoordinate3 *iterCoord3 = new SoCoordinate3;
//           iterCoord3->point.setValues(0,NUM_ITERS,iterCoords);
//           SoMarkerSet *iterMarkerSet = new SoMarkerSet;
//           iterMarkerSet->markerIndex.setValues(0,NUM_ITERS,iterMarkerIndex);
//           iterMarkerSet->markerScale.setValue( returnMarkerScale.getValue() );
//           SoMaterialBinding *iterMatBind = new SoMaterialBinding;
//           iterMatBind->value = SoMaterialBinding::PER_VERTEX;
//           SoMaterial *iterMarkerMat = new SoMaterial;
//           iterMarkerMat->diffuseColor.setValues(0,NUM_ITERS,itersColor);
//           //Add UPSTREAM to scene
//           itSep->addChild(iterCoord3);
//           itSep->addChild(iterMatBind);
//           itSep->addChild(iterMarkerMat);
//           itSep->addChild(iterMarkerSet);
//
//           //DOWNSTREAM Sub-Iterates (Marked with Open Diamonds)
//           //-------------------------------------------------------------------------------------
//           NUM_ITERS = (int) dSubIts.size();
//           SbVec3f *diterCoords = new SbVec3f[NUM_ITERS];
//           SbColor *ditersColor = new SbColor[NUM_ITERS];
//           int *diterMarkerIndex = new int[NUM_ITERS];
//           for (int i=0; i<NUM_ITERS; i++) {
//               vec42 y(0);
//               for(int j=0;j<6;j++) y[j] = dSubIts[i][j];
//               //This may need adjust for different maps
//               nvis::vec2 x = theSection.project(y).first;
//               SbVec3f mapPoint(x[0],x[1],0.0);
//               diterCoords[i].setValue(mapPoint[0],mapPoint[1],mapPoint[2]);
//               //Setting the color  - DOWNSTREAM is always the main color
//               ditersColor[i].setValue((*wCol)[0],(*wCol)[1],(*wCol)[2]);
//               diterMarkerIndex[i] = SoMarkerSet::DIAMOND_LINE_5_5;
//           }
//           SoCoordinate3 *diterCoord3 = new SoCoordinate3;
//           diterCoord3->point.setValues(0,NUM_ITERS,diterCoords);
//           SoMarkerSet *diterMarkerSet = new SoMarkerSet;
//           diterMarkerSet->markerIndex.setValues(0,NUM_ITERS,diterMarkerIndex);
//           diterMarkerSet->markerScale.setValue( returnMarkerScale.getValue() );
//           SoMaterialBinding *diterMatBind = new SoMaterialBinding;
//           diterMatBind->value = SoMaterialBinding::PER_VERTEX;
//           SoMaterial *diterMarkerMat = new SoMaterial;
//           diterMarkerMat->diffuseColor.setValues(0,NUM_ITERS,ditersColor);
//           //Add UPSTREAM to scene
//           itSep->addChild(diterCoord3);
//           itSep->addChild(diterMatBind);
//           itSep->addChild(diterMarkerMat);
//           itSep->addChild(diterMarkerSet);
//         }
//
//         //Add to a scene
//         trajScene->addChild(itSep);
//
//     }
//
//     //Show the scene
//     showGeom(trajScene);
//
//
// }
//

/// Render the manifolds
void DisplayManifolds::renderManifolds() 
{
    //Check for input data
    ManifoldData *mData = (ManifoldData*) portData.source();
    if(!mData) {
      theMsg->printf("%s: Cannot render without valid manifold data",__FILE__);
      return;
    }
    
    //Check for Connected FixedPointData
    FixedPointData *fpData = (FixedPointData*) mData->fpxDataConnection.source();
    
    if(!fpData) {
      theMsg->printf("%s: Invalid FixedPointData connected to input!",__FILE__);
      return;
    }
  
    //Clear the current scene
    scene->removeAllChildren();
    //Also clear the trajectory scene to prevent invalid references
    trajScene->removeAllChildren();

    //Map Projection:  Defines a Transform Node that applies to all rendered objects
    SoTransform* theXForm = new SoTransform();
    mapSpaceProjector.compute(horzMenu,vertMenu,horzScale,vertScale,theXForm);
    scene->addChild(theXForm);

    
    //Create a rendering queue based on display/filter options
    //-----------------------------------------------------------------------------------------
    std::queue<int> midQueueToRender;
    int numManifolds = (int) mData->theManifoldData.mapManifolds.size();
    if (filterOptions.getValue(0)==1) {
      //Only one orbit, check if we are doing only one fixed point or not
      int orbitID = filterOrbitID.getValue();
      if (filterOptions.getValue(1)==1) {
        //Only one fixed point of orbit chain
        int fpID = filterFixedPoint.getValue();
        //Select which manifolds from this orbit and particular fixed point
        switch (filterType.getValue()) {
          case 0 :
            //Only Stable manifolds
            for (int mID=0;mID<numManifolds;mID++) {
              int baseOrbitID = mData->theManifoldData.mapManifolds[mID].fpdOrbitIdx;
              int fpIndex = mData->theManifoldData.mapManifolds[mID].fpdPointIdx;
              bool isUnstable = mData->theManifoldData.mapManifolds[mID].isForward();
              if (!isUnstable && (baseOrbitID==orbitID) && (fpIndex==fpID)) 
                  midQueueToRender.push(mID);
            }
            break;
          case 1 :
            //Only Unstable manifolds
            for (int mID=0;mID<numManifolds;mID++) {
              int baseOrbitID = mData->theManifoldData.mapManifolds[mID].fpdOrbitIdx;
              int fpIndex = mData->theManifoldData.mapManifolds[mID].fpdPointIdx;
              bool isUnstable = mData->theManifoldData.mapManifolds[mID].isForward();
              if (isUnstable && (baseOrbitID==orbitID) && (fpIndex==fpID)) 
                  midQueueToRender.push(mID);
            }
            break;
          default :
            //All manifolds from orbit
            for (int mID=0;mID<numManifolds;mID++) {
              int baseOrbitID = mData->theManifoldData.mapManifolds[mID].fpdOrbitIdx;
              int fpIndex = mData->theManifoldData.mapManifolds[mID].fpdPointIdx;
              if((baseOrbitID == orbitID) && (fpIndex==fpID)) midQueueToRender.push(mID);
            }
            break;
        }
      } else {
        //Select which manifolds from this orbitID
        switch (filterType.getValue()) {
          case 0 :
            //Only Stable manifolds
            for (int mID=0;mID<numManifolds;mID++) {
              int baseOrbitID = mData->theManifoldData.mapManifolds[mID].fpdOrbitIdx;
              bool isUnstable = mData->theManifoldData.mapManifolds[mID].isForward();
              if (!isUnstable && (baseOrbitID==orbitID)) midQueueToRender.push(mID);
            }
            break;
          case 1 :
            //Only Unstable manifolds
            for (int mID=0;mID<numManifolds;mID++) {
              int baseOrbitID = mData->theManifoldData.mapManifolds[mID].fpdOrbitIdx;
              bool isUnstable = mData->theManifoldData.mapManifolds[mID].isForward();
              if (isUnstable && (baseOrbitID==orbitID)) midQueueToRender.push(mID);
            }
            break;
          default :
            //All manifolds from orbit
            for (int mID=0;mID<numManifolds;mID++) {
              int baseOrbitID = mData->theManifoldData.mapManifolds[mID].fpdOrbitIdx;
              if(baseOrbitID == orbitID) midQueueToRender.push(mID);
            }
            break;
        }
      }
        
    } else {
      //Gather manifolds based on type
      switch (filterType.getValue()) {
        case 0 :
          //Only Stable manifolds
          for (int mID=0;mID<numManifolds;mID++) {
            bool isUnstable = mData->theManifoldData.mapManifolds[mID].isForward();
            if (!isUnstable) midQueueToRender.push(mID);
          }
          break;
        case 1 :
          //Only Unstable manifolds
          for (int mID=0;mID<numManifolds;mID++) {
            bool isUnstable = mData->theManifoldData.mapManifolds[mID].isForward();
            if (isUnstable) midQueueToRender.push(mID);
          }
          break;
        default :
          //All manifolds
          for (int mID=0;mID<numManifolds;mID++) midQueueToRender.push(mID);
          break;
      }
      
    }
    
    
    
    //Render eigenVectors as line segments => NOTE: Renders ALL orbits/fixedpoints for now.
    //-----------------------------------------------------------------------------------------
    if (displayOptions.getValue(0) == 1) {
        //Unstable Eigen Vectors
        SoSeparator* eigU = new SoSeparator;
        SoLineSet* eigULine = new SoLineSet;
        SoCoordinate3* eigUCoord = new SoCoordinate3;
        SoMaterial* eigUMat = new SoMaterial;
        eigUMat->diffuseColor.setValue(1,0,1); //magenta
        //Stable Eigen Vectors
        SoSeparator* eigS = new SoSeparator;
        SoLineSet* eigSLine = new SoLineSet;
        SoCoordinate3* eigSCoord = new SoCoordinate3;
        SoMaterial* eigSMat = new SoMaterial;
        eigSMat->diffuseColor.setValue(0,1,1); //cyan

        //DrawStyle
        SoDrawStyle* eigStyle = new SoDrawStyle;
        eigStyle->lineWidth = portLineWidth.getValue();
        unsigned short pattern = 1000; //between 0 and 65535
        eigStyle->linePattern.setValue(pattern); //Dashed line

        //Add objects to scene
        eigU->addChild(eigUMat);
        eigU->addChild(eigStyle);
        eigU->addChild(eigUCoord);
        eigU->addChild(eigULine);
        scene->addChild(eigU);
        eigS->addChild(eigSMat);
        eigS->addChild(eigStyle);
        eigS->addChild(eigSCoord);
        eigS->addChild(eigSLine);
        scene->addChild(eigS);

        //Loop through all available fixed points in FixedPointData to define points
        int numfp =  fpData->getNumFixedPoints();
        int numOrbits = fpData->getNumOrbits();
        theMsg->printf("Creating eigenvectors for %d fixed points for %d orbits",numfp,numOrbits);
        if (numfp > 0) { //Don't do anything if nothing is available
            int eigUPoints = 0, eigSPoints = 0;
            int eigULineNum = 0, eigSLineNum = 0;
            //For each orbit
            for (int k=0; k<numOrbits; k++) {
              //Only look at unstable orbits (saddles)
              if (!(fpData->isSaddle(k))) continue;
              
              for (int i=0; i < fpData->getOrbitPeriod(k); i++) {
                  //First chain only here - just one orbit
                  const fixpoint& fp = fpData->getFixedPoint(k,i);
                  theMsg->printf(
                  "%s: fp %d : x0 = [%g,%g] eigVals = [%g,%g], v_Stable = [%g,%g] v_Unstable = [%g,%g]",
                  __FILE__, i, fp.pos[0],fp.pos[1],
                  fp.eval[0], fp.eval[1], fp.evec[0][0], fp.evec[0][1],
                  fp.evec[1][0], fp.evec[1][1]);
                  //Add Unstable coords & line
                  double s = (double) eigenLength.getValue();
                  eigUCoord->point.set1Value(eigUPoints++,
                          fp.pos[0]-s*fp.evec[1][0],
                          fp.pos[1]-s*fp.evec[1][1],
                          0.0); //-scale*eigVec
                  eigUCoord->point.set1Value(eigUPoints++,
                          fp.pos[0]+s*fp.evec[1][0],
                          fp.pos[1]+s*fp.evec[1][1],
                          0.0); //+scale*eigVec
                  eigULine->numVertices.set1Value(eigULineNum,2); //Add the two-point line
                  eigULineNum++;
                  //Add Stable Coords & line
                  eigSCoord->point.set1Value(eigSPoints++,
                          fp.pos[0]-s*fp.evec[0][0],
                          fp.pos[1]-s*fp.evec[0][1],
                          0.0); //-scale*eigVec
                  eigSCoord->point.set1Value(eigSPoints++,
                          fp.pos[0]+s*fp.evec[0][0],
                          fp.pos[1]+s*fp.evec[0][1],
                          0.0); //+scale*eigVec
                  eigSLine->numVertices.set1Value(eigSLineNum,2); //Add the two-point line
                  eigSLineNum++;
              }
           }
        }
    } //End Eigenvector Render

    
    //-----------------------------------------------------------------------------------------
    //Render manifolds
    //-----------------------------------------------------------------------------------------
    SoSeparator* manifoldSep = new SoSeparator;
    SoSeparator* mSep = new SoSeparator; //Unstable
    SoLineSet* mLines = new SoLineSet;
    SoCoordinate3* mCoords = new SoCoordinate3;
    SoSeparator* smSep = new SoSeparator; //Stable
    SoLineSet* smLines = new SoLineSet;
    SoCoordinate3* smCoords = new SoCoordinate3;
    manifoldSep->addChild(mSep);
    manifoldSep->addChild(smSep);
    targetNode = manifoldSep;
    wuTarget = mSep;
    wsTarget = smSep;
    scene->addChild(manifoldSep);

    //SoDrawStyle -> change the way lines are rendered
    SoDrawStyle* manifoldStyle = new SoDrawStyle;
    manifoldStyle->lineWidth = portLineWidth.getValue();

    //Materials for manifold lines
    SoMaterialBinding *mMatBind = new SoMaterialBinding;
    //mMatBind->value = SoMaterialBinding::PER_PART;
    mMatBind->value = SoMaterialBinding::PER_VERTEX; //per-vertex color
    SoMaterial *mMat = new SoMaterial; //Unstable
    McColor *mCol = new McColor;
    colors.getColor(*mCol,2);
    mMat->diffuseColor.setValue((*mCol)[0],(*mCol)[1],(*mCol)[2]);
    SoMaterialBinding *smMatBind = new SoMaterialBinding;
    //smMatBind->value = SoMaterialBinding::PER_PART; //A color for each Segment
    smMatBind->value = SoMaterialBinding::PER_VERTEX; //Vertexes
    SoMaterial *smMat = new SoMaterial; //Stable
    McColor *smCol = new McColor;
    colors.getColor(*smCol,1);
    smMat->diffuseColor.setValue((*smCol)[0],(*smCol)[1],(*smCol)[2]);

    //Determine if we will show the manifold objects:
    if(displayOptions.getValue(3) == 1) {
      //Add objects to scene graph
      mSep->addChild(mMatBind);
      mSep->addChild(mMat);
      mSep->addChild(manifoldStyle);
      mSep->addChild(mCoords);
      mSep->addChild(mLines);
      smSep->addChild(smMatBind);
      smSep->addChild(smMat);
      smSep->addChild(manifoldStyle);
      smSep->addChild(smCoords);
      smSep->addChild(smLines);
    }

    //Fill Coordinate Arrays for each Manifold object (from rendering queue)
    typedef std::pair<int,int> IntPair;
    uLinesetToManifold.clear(); sLinesetToManifold.clear();
    int unstableSegTotal = 0, unstablePoints = 0;
    int stableSegTotal = 0, stablePoints = 0;
    int uLineNum = 0;
    int sLineNum = 0;
    //Make the main render queue as a priority_queue that puts long, highly unstable orbits first
    // (putting them on the bottom of the visualization).
    ManifoldDataRenderPriority renderPriority(mData->theManifoldData);
    std::priority_queue<int,std::vector<int>,ManifoldDataRenderPriority> mainRenderQueue(renderPriority);
    //Emplace queue objects
    while(!midQueueToRender.empty()) {
      mainRenderQueue.push(midQueueToRender.front());
      midQueueToRender.pop();
    }
    //For each manifold object in queue
    while (!mainRenderQueue.empty()) {
        int mID = mainRenderQueue.top(); //Priority rendering
        int ptsPerLine = 0;
        bool isUnstable = mData->theManifoldData.mapManifolds[mID].isForward();
        //theMsg->printf("Starting new %s manifold...",(isUnstable)? "UNSTABLE" : "STABLE");
        nvis::vec2 prev(0.0,0.0);
        bool firstSeg = true, connected = true;
        //int fOrbitID = fpData->getOrbitID(mData->theManifoldData.mapManifolds[mID].fpdOrbitIdx,
        //                                  mData->theManifoldData.mapManifolds[mID].fpdPointIdx);
        int fOrbitID = mData->theManifoldData.mapManifolds[mID].fpdOrbitIdx;
        //For each segment
        std::vector<ManifoldSeg>::iterator segit;
        for ( segit =  mData->theManifoldData.mapManifolds[mID].segments.begin();
              segit != mData->theManifoldData.mapManifolds[mID].segments.end(); segit++) {
              //The two points of this segment
              nvis::vec2& x0 = (*segit)[0];
              nvis::vec2& x1 = (*segit)[1];
              IntPair manIntPair(mID,segit->segID);
        
              //Check if this point will be within bounds
              bool inBounds = true;
              nvis::bbox2 bbox(nvis::vec2(portXbounds.getValue(0),portYbounds.getValue(0)),
                               nvis::vec2(portXbounds.getValue(1),portYbounds.getValue(1)) );
              if ( !bbox.inside(x0) || !bbox.inside(x1) ) {
                connected = false;
                inBounds = false;
              }
              if (!inBounds) continue; //Go to next segment
        
              //Is this connected to the last segment
              if (!firstSeg) {
                 if(nvis::all( x0 == prev )) {
                   connected = true;
                 } else {
                   connected = false;
                 }
              }
              //If first segment or not connected, add both points
              if(firstSeg || !connected) {
                //If not connected, we must end the last line
                if(!connected) {
                  if(isUnstable) {
                    mLines->numVertices.set1Value(uLineNum,ptsPerLine);
                    uLineNum++;
                  } else {
                    smLines->numVertices.set1Value(sLineNum,ptsPerLine);
                    sLineNum++;
                  }
                }
                if(isUnstable) {
                  //Add points/color to lineset containers
                  mCoords->point.set1Value(unstablePoints++,x0[0],x0[1],0.0);
                  mMat->diffuseColor.set1Value(unstablePoints,unstableColormap.getColor(fOrbitID));
                  mCoords->point.set1Value(unstablePoints++,x1[0],x1[1],0.0);
                  mMat->diffuseColor.set1Value(unstablePoints,unstableColormap.getColor(fOrbitID));
                  //theMsg->printf("Debug:  Adding NEW UNSTABLE segment %d : [%g %g] - [%g %g]",unstableSegTotal,x0[0],x0[1],x1[0],x1[1]);
                  //Add to index map
                  LineSegPair thisSeg(uLineNum,unstableSegTotal);
                  uLinesetToManifold.insert( std::pair<LineSegPair,IntPair>(thisSeg,manIntPair) );
                  unstableSegTotal++;
                } else {
                  //Add points/color to lineset containers
                  smCoords->point.set1Value(stablePoints++,x0[0],x0[1],0.0);
                  smMat->diffuseColor.set1Value(stablePoints,stableColormap.getColor(fOrbitID));
                  smCoords->point.set1Value(stablePoints++,x1[0],x1[1],0.0);
                  smMat->diffuseColor.set1Value(stablePoints,stableColormap.getColor(fOrbitID));
                  //theMsg->printf("Debug:  Adding NEW STABLE segment %d : [%g %g] - [%g %g]",stableSegTotal,x0[0],x0[1],x1[0],x1[1]);
                  //Add to index map
                  LineSegPair thisSeg(sLineNum,stableSegTotal);
                  sLinesetToManifold.insert( std::pair<LineSegPair,IntPair>(thisSeg,manIntPair) );
                  stableSegTotal++;                  
                }
                
                //Update to no longer be first segment
                if(firstSeg) firstSeg = false;
                //Restart the pointsPerLine Counter
                ptsPerLine = 2;
              } else { //Already connected so just add a single point
                if(isUnstable) {
                  mCoords->point.set1Value(unstablePoints++,x1[0],x1[1],0.0);
                  mMat->diffuseColor.set1Value(unstablePoints,unstableColormap.getColor(fOrbitID));
                  //theMsg->printf("Debug:  Adding UNSTABLE segment %d : [%g %g] - [%g %g]",unstableSegTotal,x0[0],x0[1],x1[0],x1[1]);
                  //Add to index map
                  LineSegPair thisSeg(uLineNum,unstableSegTotal);
                  uLinesetToManifold.insert( std::pair<LineSegPair,IntPair>(thisSeg,manIntPair) );
                  unstableSegTotal++;
                } else {
                  smCoords->point.set1Value(stablePoints++,x1[0],x1[1],0.0);
                  smMat->diffuseColor.set1Value(stablePoints,stableColormap.getColor(fOrbitID));    
                  //theMsg->printf("Debug:  Adding NEW STABLE segment %d : [%g %g] - [%g %g]",stableSegTotal,x0[0],x0[1],x1[0],x1[1]);
                  //Add to index map
                  LineSegPair thisSeg(sLineNum,stableSegTotal);
                  sLinesetToManifold.insert( std::pair<LineSegPair,IntPair>(thisSeg,manIntPair) );
                  stableSegTotal++;    
                }
                ptsPerLine++;
              }
              //Store the last point to see if we need it again
              prev = x1;
        }
        //Always close the last line
        if(isUnstable) {
            mLines->numVertices.set1Value(uLineNum,ptsPerLine);
            uLineNum++;
        } else {
            smLines->numVertices.set1Value(sLineNum,ptsPerLine);
            sLineNum++;
        }
        
        //Progress queue
        mainRenderQueue.pop();
        
    } //End while loop for main render
    
    //Output
    /*theMsg->printf("%s: Stable Lines = %d with %d segments",__FILE__,sLineNum,stableSegTotal);
    theMsg->printf("%s: Elements in StableLineSet Map = %d",__FILE__,(int)sLinesetToManifold.size());
    theMsg->printf("%s: Unstable Lines = %d with %d segments",__FILE__,uLineNum,unstableSegTotal);
    theMsg->printf("%s: Elements in UnstableLineSet Map = %d",__FILE__,(int)uLinesetToManifold.size());
    */


    //Render fixed points - Only currently renders ALL
    //---------------------------------------------------------------------------------------
    if (displayOptions.getValue(1) == 1) { 
        //After manifolds to make them appear "on top"
        SoSeparator *saddleSep = new SoSeparator;
        SoPointSet *sSet = new SoPointSet;
        SoDrawStyle *drawStyle = new SoDrawStyle;
        drawStyle->pointSize = portPointSize.getValue();

        //Materials
        SoVertexProperty *sProp = new SoVertexProperty;
        McColor *sCol = new McColor;
        colors.getColor(*sCol,0);
        sProp->orderedRGBA = colors.getPackedColor(0);


        //Fill Coordinate arrays
        int nSaddles = -1;
        int numOrbits = fpData->getNumOrbits();
        std::vector<SbVec3f> verts;
        for(int k=0;k<numOrbits;k++) {
          if(!fpData->isSaddle(k)) continue;
          
          for(int i=0;i<fpData->getOrbitPeriod(k); i++) {
            const xavier::fixpoint& fp = fpData->getFixedPoint(k,i);
            SbVec3f vertexPos;
            vertexPos[0] = fp.pos[0];
            vertexPos[1] = fp.pos[1];
            vertexPos[2] = 0.0;
            verts.push_back( vertexPos );
          }
        }
        nSaddles = (int)verts.size();
        SbVec3f *vertexPos = new SbVec3f[nSaddles];
        for(int i=0;i<nSaddles;i++) vertexPos[i] = verts[i];
        sProp->vertex.setValues(0,nSaddles,vertexPos);
        //Add to scene
        sSet->vertexProperty.setValue(sProp);
        saddleSep->addChild(drawStyle);
        saddleSep->addChild(sSet);
        saddleSep->addChild(sProp);
        scene->addChild(saddleSep);
    }


    //---------------------------------------------------------------------------------------
    //Render the Map Discontinuities
    //---------------------------------------------------------------------------------------
    if (displayOptions.getValue(2)==1) {
      SoSeparator *mapDisSep = new SoSeparator;
      SoDepthOffset *mapDisDOff = new SoDepthOffset;
      SoCoordinate3 *mapDisCoords = new SoCoordinate3;
      SoDrawStyle *mapDisDraw = new SoDrawStyle;
      mapDisDraw->pointSize = portPointSize.getValue();
      SoMarkerSet *mapDisMarkers = new SoMarkerSet;
      SoMaterialBinding *mapDisBinder = new SoMaterialBinding;
      mapDisBinder->value = SoMaterialBinding::PER_VERTEX;
      SoMaterial *mapDisMat = new SoMaterial;

      //Iterate through the sepList and grab all points from manifold
      std::vector<nvis::vec2> xFail;
      std::vector<SbColor> failColors;
      std::list< MapDiscont >::const_iterator failIT;
      int failPtIdx = 0;
      std::queue<int> renderQueue(midQueueToRender);
      //For each Manifold in render queue
      while (!renderQueue.empty()) {
        int k = renderQueue.front();
        std::list<MapDiscont> sepList; 
        mData->theManifoldData.mapManifolds[k].getMapDiscontinuityList(sepList);
        //For each MapDiscont
        for(failIT = sepList.begin(); failIT!=sepList.end(); ++failIT) {
          xFail.push_back( failIT->where() );
          McColor failColor;
          switch ( failIT->type ) {
            case MapDiscont::BACKWARD_MAP :
              mapDisMarkers->markerIndex.set1Value(failPtIdx,SoMarkerSet::PINE_TREE_FILLED_9_9);
              failColor.setValue(0.4,0.4,1);
              break;
            case MapDiscont::FORWARD_MAP :
              mapDisMarkers->markerIndex.set1Value(failPtIdx,SoMarkerSet::SHIP_FILLED_9_9);
              failColor.setValue(1,0.4,0.4);
              break;
            case MapDiscont::BACKWARD_SINGULARITY :
              mapDisMarkers->markerIndex.set1Value(failPtIdx,SoMarkerSet::CIRCLE_LINE_9_9);
              failColor.setValue(0,1,1);
              break;
            case MapDiscont::FORWARD_SINGULARITY :
              mapDisMarkers->markerIndex.set1Value(failPtIdx,SoMarkerSet::CROSS_9_9);
              failColor.setValue(1,0,1);
              break;
            case MapDiscont::BACKWARD_SECTION_SEP :
              //Blue rhombus (filled)
              mapDisMarkers->markerIndex.set1Value(failPtIdx,SoMarkerSet::RHOMBUS_LINE_9_9);
              failColor.setValue(0,0,1);
              break;
            case MapDiscont::FORWARD_SECTION_SEP :
              //Red rhombus
              mapDisMarkers->markerIndex.set1Value(failPtIdx,SoMarkerSet::RHOMBUS_LINE_9_9);
              failColor.setValue(1,0,0);
              break;
            case MapDiscont::BACKWARD_SECTION_SEP_NODE :
              //Blue rhombus (filled)
              mapDisMarkers->markerIndex.set1Value(failPtIdx,SoMarkerSet::DIAMOND_LINE_9_9);
              failColor.setValue(0,0,1);
              break;
            case MapDiscont::FORWARD_SECTION_SEP_NODE :
              //Red rhombus
              mapDisMarkers->markerIndex.set1Value(failPtIdx,SoMarkerSet::DIAMOND_LINE_9_9);
              failColor.setValue(1,0,0);
              break;
            case MapDiscont::DOUBLE_PERIOD_FIXED_POINT :
              //Gray Y
              mapDisMarkers->markerIndex.set1Value(failPtIdx,SoMarkerSet::Y_9_9);
              failColor.setValue(0.6,0.6,0.6);
              break;
            case MapDiscont::FIXED_POINT_SUSPECTED :
              //Gray star
              mapDisMarkers->markerIndex.set1Value(failPtIdx,SoMarkerSet::STAR_9_9);
              failColor.setValue(0.6,0.6,0.6);
              break;
            case MapDiscont::UNRESOLVED :
              //Black Hourglass
              mapDisMarkers->markerIndex.set1Value(failPtIdx,SoMarkerSet::HOURGLASS_FILLED_9_9);
              if (failIT->period < 0) {
                failColor.setValue(0.2,0.2,0.4);
              } else
                failColor.setValue(0.4,0.2,0.2);
              break;
            default :
              //Unknown failure
              mapDisMarkers->markerIndex.set1Value(failPtIdx,SoMarkerSet::CAUTION_LINE_9_9);
              failColor.setValue(1,1.0,0.0);
              break;
          }
          //If fixed, force marker color to green
          //if (failIT->fixed) failColor.setValue(0,1,0);
          //Add Marker Colors
          failColors.push_back(SbColor(failColor.r,failColor.g,failColor.b));
          //failColors.push_back(failColor.getPackedColor(0));
          //Increment counter
          failPtIdx++;

        } //End SepList Loop
        
        //Mark the last seed point - End of processed manifold
        nvis::vec2 lastSeedPt = (mData->theManifoldData.mapManifolds[k].getWorkingSegment())[1];
        xFail.push_back( lastSeedPt );
        mapDisMarkers->markerIndex.set1Value(failPtIdx,SoMarkerSet::BAR_9_9);
        failColors.push_back(SbColor(0,1,0));
        failPtIdx++;

        //Move to next manifold
        renderQueue.pop();
      } //End While Loop
      
      //Collect all data for rendering
      int numFailPoints = (int) failColors.size();
      theMsg->printf("There are %d map discontinuities from this manifold",numFailPoints);
      SbColor *failPointColors = new SbColor[numFailPoints];
      SbVec3f *failPts = new SbVec3f[numFailPoints];
      for(int i=0;i<numFailPoints;i++) {
        failPts[i] = SbVec3f(xFail[i][0],xFail[i][1],0.0);
        failPointColors[i] = failColors[i];
      }
      mapDisCoords->point.setValues(0,numFailPoints,failPts);
      mapDisMat->diffuseColor.setValues(0,numFailPoints,failPointColors);

      //Add to scene
      mapDisSep->addChild(mapDisDOff);
      mapDisSep->addChild(mapDisDraw);
      mapDisSep->addChild(mapDisCoords);
      mapDisSep->addChild(mapDisBinder);
      mapDisSep->addChild(mapDisMat);
      mapDisSep->addChild(mapDisMarkers);
      scene->addChild(mapDisSep);

    }
    

    //Insert the pick-callback FIRST into the local scene (for interaction)
    //scene->insertChild(eventCB,0);
    //Actually need at root-level, have to work with multiple scenes
    
    //Render the scene
    showGeom(scene);
}

/// Display the topology design scene if input IState available (constraining line with dV values)
void DisplayManifolds::renderTopoDesignScene()
{
    //Make this scene full of small objects as it is continuously updated
    tdScene->removeAllChildren();
    //Render on top of stuff
    SoDepthOffset *tdOffset = new SoDepthOffset();
    tdScene->addChild(tdOffset);
    
    //Manifold and Fixed point data (will be avaliable)
    ManifoldData *mData = (ManifoldData*) portData.source();
    FixedPointData *fpData = (FixedPointData*) mData->fpxDataConnection.source();
    
    //Input state : mark map location (x,xdot)
    HxInitialState* iState = (HxInitialState*) initStateConnection.source();
    //Could Input be invalid here?
    if (!iState) return;
    
    //Initial state data:
    State &theState = iState->data.ic;
    SbVec3f iStateMapCoord(theState.x,theState.xd,0.0); //Single precision
    //theMsg->printf("RenderTDS: istate Map = (%g, %g)", theState.x, theState.xd);
    double inputC = sys->getJacobiConstant(theState);
    
    //Map position for picked point OR selected manifold
    nvis::vec2 xPP(0.0);
    
    //Place-holder value for orbit ID
    int fpOrbitID = 0;
    
    //ONLY Render Rest if we have a valid design selection
    if (topoDesignSelected) {
        ///Convert segment, line index, and parameter into manifold and state
        LineSegPair pickedManifold( selectedLineID, selectedPartID );
        //LineSegPair hcPartnerManifold( hcPartnerLineID, hcPartnerPartID );
        //LineID,SegID -> ManifoldID
        IntPair manSegPair;
        //IntPair hcManSegPair;
        if(selectedUnstable) {
          //Get (manifoldID,segmentID)
          manSegPair = uLinesetToManifold[pickedManifold];
        } else {
          manSegPair = sLinesetToManifold[pickedManifold];
        }
        /*//Heteroclinic piece - Ignore
        if(hcOutput.getValue(0)==1) {
          if(selectedUnstable) {
            //Get (manifoldID,segmentID) for HC partner
            hcManSegPair = sLinesetToManifold[hcPartnerManifold];
          } else {
            hcManSegPair = uLinesetToManifold[hcPartnerManifold];
          }
        }*/
        
        //Get the Manifold and Segment 
        int mID = manSegPair.first; int segID = manSegPair.second;
        const ManifoldData::ManifoldType& theManifold = mData->theManifoldData.mapManifolds[mID];
        const ManifoldData::ManifoldSeg& theSeg = theManifold.segments[segID];
        fpOrbitID = theManifold.fpdOrbitIdx;
        
        //Construct a dummy manifold segment to hold the deltaV-line information
        ManifoldData::ManifoldSeg dVLineSeg( 
            nvis::vec2(theState.x,(double)portYbounds.getValue(0)),
            nvis::vec2(theState.x,(double)portYbounds.getValue(1)),
            nvis::vec3(0.0) ); //Dummy data [not used]
        //Compute the linear parameter for selected point
        //selectedLinearParameter = theSeg.getLinearParam(xApprox);
        
        //Skip because already computed in topologyDeltaV()
        //Compute parameter based on interesection of deltaV line and the selected segment
        /*if ( !theSeg.intersect(dVLineSeg) ) {
          theMsg->printf("%s : Warning! DeltaV Line and selected manifold do not intersect!",getName());
        }
        double dummyTau = 0.0;
        theSeg.getIntersection(dVLineSeg,selectedLinearParameter,dummyTau);*/
        //The actual selected manifold point to double precision:
        xPP = theSeg.getPoint(selectedLinearParameter); //MapState (x,xdot)
        
    }
    else {
      //Otherwise, just use the constrainted picked point to compute the deltaV
      SbVec3f mp;
      mp = mapSpaceProjector.worldToMap(horzMenu,vertMenu,horzScale,vertScale,pickedPoint);
      xPP[0] = mp[0]; xPP[1] = mp[1]; //Current pickedPoint location on map 
    }
    
    
    //CR3BP info
    double mup = sys->mup;
    double jcManifolds = fpData->getJacobiConstant();
    rhs_type   theEOMs(jcManifolds,mup);  //The right-hand side:  CR3BP EoMs
    section_type  theSection(theEOMs); //Hyperplane: y=0 by planar_section

    // Get the full state of the manifold with a reverse projection call
    State ppState;
    vec42 y = theSection.unproject(xPP);
    for(int i=0;i<6;i++) ppState[i] = y[i];
    
    //Transform 
    SoTransform *mapToWorldTrans = new SoTransform;
    mapSpaceProjector.getTransform(mapToWorldTrans);
    tdScene->addChild(mapToWorldTrans);
    
    //Picked Point:
    //SbVec3f mapPoint = mapSpaceProjector.worldToMap(
    //                     horzMenu,vertMenu,horzScale,vertScale,
    //                     pickedPoint);
    SbVec3f mapPoint(xPP[0],xPP[1],0.0); //Single precision

        
    //Compute the DeltaV (Before rendering to check what's valid)
    bool validPoint = true;
    double dxd = 0.0, dyd = 0.0, dvND = 0.0, dC = 0.0, disc = 0.0;
    bool fwd = true;
    //Split settings based on design selection :
    //--------------------------------------------------------------------------------------
    //State to Ws (F) || State to Wu (B)
    if (portDesignOpt.getValue() == 0 || portDesignOpt.getValue() == 2) {
        //Set reverse design direction if going backward
        if(portDesignOpt.getValue() == 2) fwd = false;
        //Translation 
        dxd = xPP[1] - theState.xd;
        //Check sign
        if (!fwd) dxd *= -1.0;
        //Mark C0 value for initial state 
        dC = jcManifolds - inputC;
        //MapManeuver: InitialState to picked point
        disc = theState.yd*theState.yd - 2.0*theState.xd*dxd-dxd*dxd;
        //Only apply dC change when large enough
        if (std::fabs(dC)>= 1e-8) {
            //Free form definition
            disc = theState.yd*theState.yd - 2.0*theState.xd*dxd-dxd*dxd - dC;
        } 
        //Check zero-velocity condition, only proceed if valid 
        if (disc < 0.0) {
            //Set the selection point as red, indicating invalid!
            validPoint = false;
            //Keep values at zero.
        } else {
            //Compute the cheapest MAGNITUDE
            double dydot0 = -theState.yd + sqrt(disc);
            double dydot1 = -theState.yd - sqrt(disc);
            dyd = (std::fabs(dydot0)<std::fabs(dydot1))? dydot0 : dydot1;
            dvND = sqrt(dxd*dxd + dyd*dyd);
        }
    }     
    //--------------------------------------------------------------------------------------
    //Wu to State (F) || Ws to State (B)
    else if (portDesignOpt.getValue() == 1 || portDesignOpt.getValue() == 3) { 
        //Set reverse design direction if going backward
        if(portDesignOpt.getValue() == 3) fwd = false;
        //Translation
        dxd = theState.xd - xPP[1];
        //Check sign!!!
        if (!fwd) dxd *= -1.0;
        //Jacobi change is from selection to input state  
        dC = inputC - jcManifolds;
        //MapManeuver: InitialState to picked point
        disc = ppState.yd*ppState.yd - 2.0*ppState.xd*dxd-dxd*dxd;
        //Only apply dC change when large enough
        if (std::fabs(dC)>= 1e-8) {
            //Free form definition
            disc = ppState.yd*ppState.yd - 2.0*ppState.xd*dxd-dxd*dxd - dC;
        } 
        //Check zero-velocity condition, only proceed if valid 
        if (disc < 0.0) {
            //Set the selection point as red, indicating invalid!
            validPoint = false;
            //Keep values at zero.
        } else {
            //Compute the cheapest MAGNITUDE
            double dydot0 = -ppState.yd + sqrt(disc);
            double dydot1 = -ppState.yd - sqrt(disc);
            dyd = (std::fabs(dydot0)<std::fabs(dydot1))? dydot0 : dydot1;
            dvND = sqrt(dxd*dxd + dyd*dyd);
        }    
    }
    //theMsg->printf("RenderTDS: DeltaV Info:  %s", ((validPoint)?"Valid":"INVALID"));
    //theMsg->printf("  dC = %g  dvND = %g  with [%g %g %g]",dC,dvND,dxd,dyd,0.0);
    //theMsg->printf("  iStateMapCoord = (%g, %g, %g)", 
    //       iStateMapCoord[0],iStateMapCoord[1],iStateMapCoord[2]);
    
    //---------------------------------------------------------------------------------------
    // Render the Topology Design Line constraint (Vertical translation)
    //---------------------------------------------------------------------------------------
    //Build a line to confine the translation of a selection (vertical movement in xdot)
    //Use bounds to figure out how long to make line
    float constraintVerts[2][3] = {
        {iStateMapCoord[0], portYbounds.getValue(0), 0.0},
        {iStateMapCoord[0], portYbounds.getValue(1), 0.0}
    };  
    SoSeparator *lsSep = new SoSeparator;
    SoCoordinate3 *lsCoords = new SoCoordinate3;
    SoLineSet *lsLine = new SoLineSet;
    //Draw Style
    SoDrawStyle *lStyle = new SoDrawStyle;
    lStyle->lineWidth = portLineWidth.getValue();
    lStyle->linePattern.setValue( 10110 );
    SoMaterialBinding *lsMatBind = new SoMaterialBinding;
    //One color for each line segment
    lsMatBind->value = SoMaterialBinding::PER_PART;
    SoMaterial *lsMat = new SoMaterial;
    lsCoords->point.setValues(0,2,constraintVerts);
    lsLine->numVertices.set1Value(0,2);
    McColor cLineCol; designColors.getColor(cLineCol,2);
    lsMat->diffuseColor.set1Value(0,SbColor(cLineCol.r,cLineCol.g,cLineCol.b));
    //Assign design selection node
    topoDesignTargetNode = lsSep;
    //Add to display
    lsSep->addChild(lsMatBind);
    lsSep->addChild(lsMat);
    lsSep->addChild(lStyle);
    lsSep->addChild(lsCoords);
    lsSep->addChild(lsLine);
    tdScene->addChild(lsSep);
     
    
    //Display the initial state and selection point 
    SbVec3f *markerCoord = new SbVec3f[2];
    markerCoord[0] = iStateMapCoord;
    markerCoord[1] = SbVec3f(iStateMapCoord[0],mapPoint[1],0.0);
    int markerIndex[2] = { 
        SoMarkerSet::CIRCLE_FILLED_9_9,
        SoMarkerSet::SATELLITE_FILLED_9_9
    };
    SbColor *markerColor = new SbColor[2];
    McColor iStateCol; designColors.getColor(iStateCol,0);
    markerColor[0].setValue(iStateCol.r,iStateCol.g,iStateCol.b);
    McColor ppCol; designColors.getColor(ppCol,1);
    //Only color working design with correct color:
    if (validPoint) {
        markerColor[1].setValue(ppCol.r,ppCol.g,ppCol.b);
    } else {
        markerColor[1].setValue(1.0,0.0,0.0);
    }
    SoCoordinate3 *markerCoord3 = new SoCoordinate3;
    markerCoord3->point.setValues(0,2,markerCoord);
    SoMarkerSet *markerSet = new SoMarkerSet;
    markerSet->markerIndex.setValues(0,2,markerIndex);
    markerSet->markerScale.setValue( portMarkerScale.getValue() );
    SoMaterialBinding *matBind = new SoMaterialBinding;
    matBind->value = SoMaterialBinding::PER_VERTEX;
    SoMaterial *markerMat = new SoMaterial;
    markerMat->diffuseColor.setValues(0,2,markerColor);
    tdScene->addChild(markerCoord3);
    tdScene->addChild(matBind);
    tdScene->addChild(markerMat);
    tdScene->addChild(markerSet);
      
    //Make sure initial state is within bounds selected by object
      // ->Skip for now...  assume ok...
    
    //Display the translation from an input point as vector 
    float dvVerts[2][3] = {
        {iStateMapCoord[0], iStateMapCoord[1], 0.0},
        {iStateMapCoord[0], mapPoint[1], 0.0}
    };  
    SoSeparator *lsSep2 = new SoSeparator;
    SoCoordinate3 *lsCoords2 = new SoCoordinate3;
    SoLineSet *lsLine2 = new SoLineSet;    
    //Draw Style
    SoDrawStyle *lsStyle2 = new SoDrawStyle;
    lsStyle2->lineWidth = 2.0*portLineWidth.getValue();
    SoMaterialBinding *lsMatBind2 = new SoMaterialBinding;
    //One color for each line segment
    lsMatBind2->value = SoMaterialBinding::PER_PART;
    SoMaterial *lsMat2 = new SoMaterial;
    lsCoords2->point.setValues(0,2,dvVerts);
    lsLine2->numVertices.set1Value(0,2);
    McColor dvLineCol; designColors.getColor(dvLineCol,1);
    lsMat2->diffuseColor.set1Value(0,SbColor(dvLineCol.r,dvLineCol.g,dvLineCol.b));
    lsSep2->addChild(lsMatBind2);
    lsSep2->addChild(lsMat2);
    lsSep2->addChild(lsStyle2);
    lsSep2->addChild(lsCoords2);
    lsSep2->addChild(lsLine2);
    tdScene->addChild(lsSep2);

    //Render the deltaV text
    if (portDesignDisplayOpts.getValue(0) == 1) { 
      //Create a SoText2 to display the velocity magnitude
      SbVec3f dvTextMapCoord(iStateMapCoord[0],(iStateMapCoord[1]+mapPoint[1])/2.0,0.0);
      SoMaterial *ppMat = new SoMaterial;
      if(validPoint) {
          McColor fontCol( portDesignFont.getFontColor() );
          ppMat->diffuseColor.setValue(fontCol.r,fontCol.g,fontCol.b);
      } else {
          //Indicate an invalid point by red!
          ppMat->diffuseColor.setValue(SbColor(1.0,0.0,0.0));
      }
      //Font
      SoFont *textFont = new SoFont;
      QString familyName = portDesignFont.getFontName();
      bool isBold = portDesignFont.isBoldFont();
      bool isItalic = portDesignFont.isItalicFont();
      if (isBold && isItalic) {
        familyName.append(" : Bold Italic");
      } else if (isBold) {
        familyName.append(" : Bold");
      } else if (isItalic) {
        familyName.append(" : Italic");
      }
      QByteArray fbarray = familyName.toLocal8Bit();
      textFont->name.setValue(fbarray.data());
      textFont->size.setValue(portDesignFont.getFontSize());
      //Add to sep to control all subsequent SoNodes
      SoSeparator *textSep = new SoSeparator;
      //textSep->addChild(textDOffset);
      textSep->addChild(ppMat);
      textSep->addChild(textFont);
      //Indicate a translation to define the location 
      SoTransform *textXform = new SoTransform;
      textXform->translation.setValue( dvTextMapCoord - SbVec3f(0,0,0) );
      //The actual Text 
      QString dvTextStr, dxdTextStr;
      double value = 0.0, value2 = 0.0;
      if (validPoint) {
          switch(portDeltaVUnits.getValue()) {
              case 1: //km/s 
                  //Convert the velocity magnitude to km/s 
                  value = sys->getDimVelocity(dvND);
                  dvTextStr.setNum(value,'g',5);
                  dvTextStr += " km/s";
                  value2 = sys->getDimVelocity(dxd);
                  dxdTextStr.setNum(value2,'g',5);
                  dxdTextStr += " km/s";
                  break;
              case 2: //m/s
                  //Convert the velocity magnitude to km/s 
                  value = sys->getDimVelocity(dvND);
                  dvTextStr.setNum(value*1.e3,'g',5);//m/s
                  dvTextStr += " m/s";
                  value2 = sys->getDimVelocity(dxd);
                  dxdTextStr.setNum(value2*1.e3,'g',5);
                  dxdTextStr += " m/s";
                  break;
              default : //ND 
                  dvTextStr.setNum(dvND,'g',5);
                  dvTextStr += " ND";
                  dxdTextStr.setNum(dxd,'g',5);
                  dxdTextStr += " ND";
                  break;
          }
      } else {
          dvTextStr = "N/A";
      }
      //Move off the line a bit and add vertical translation:
      QString fullText;
      if (portDesignDisplayOpts.getValue(1)==1) {
        fullText = " " + QString(QChar(0x0394)) + "V= " + dvTextStr;
      }
      if (validPoint && portDesignDisplayOpts.getValue(2)==1) { 
        fullText += "  [" + QString(QChar(0x0394)) + "xd= " + dxdTextStr + "]";
      }
      SoText2 *dvText = new SoText2;      
      SbString soStr; soStr.fromUtf16(fullText.utf16());
      dvText->string.setValue(soStr);
      //Add to scene 
      textSep->addChild(textXform);
      textSep->addChild(dvText);
      tdScene->addChild(textSep);
    }
    
    //Display the OrbitID of the selection
    if (topoDesignSelected && portDesignDisplayOpts.getValue(1)==1) {
      //Create a SoText2 to display orbit ID of the selected manifold
      SbVec3f oidMapCoord(xPP[0],xPP[1],0.0);
      //Get color from colormap
      SoMaterial *oidMat = new SoMaterial;
      if(selectedUnstable) {
          oidMat->diffuseColor.setValue(unstableColormap.getColor(fpOrbitID));
      } else {
          oidMat->diffuseColor.setValue(stableColormap.getColor(fpOrbitID));
      }
      //Font
      SoFont *textFont = new SoFont;
      QString familyName = portDesignFont.getFontName();
      bool isBold = portDesignFont.isBoldFont();
      bool isItalic = portDesignFont.isItalicFont();
      if (isBold && isItalic) {
        familyName.append(" : Bold Italic");
      } else if (isBold) {
        familyName.append(" : Bold");
      } else if (isItalic) {
        familyName.append(" : Italic");
      }
      QByteArray fbarray = familyName.toLocal8Bit();
      textFont->name.setValue(fbarray.data());
      textFont->size.setValue(portDesignFont.getFontSize());
      
      //Add to sep to control all subsequent SoNodes
      SoSeparator *textSep = new SoSeparator;
      //textSep->addChild(textDOffset);
      textSep->addChild(oidMat);
      textSep->addChild(textFont);
      //Indicate a translation to define the location 
      SoTransform *textXform = new SoTransform;
      textXform->translation.setValue( oidMapCoord - SbVec3f(0,0,0) );
      //The actual Text 
      QString oidTextStr;
      oidTextStr.setNum( fpOrbitID );
      //Move off the line a bit and add vertical translation:
      QString fullText = " ID = " + oidTextStr;
      SoText2 *oidText = new SoText2;
      oidText->string.setValue(fullText.toLocal8Bit().data());
      //Add to scene 
      textSep->addChild(textXform);
      textSep->addChild(oidText);
      tdScene->addChild(textSep);
    }
    
    
    //Show the scene 
    showGeom(tdScene);
}



/// Render the translation (or DeltaV Line) during selection for WWDesign
void DisplayManifolds::renderWWTempScene(bool remove)
{
    //Reset the scene
    wwdvScene->removeAllChildren();
    if(remove) return;
    
    //Don't do this if there is an invalid design start
    if(!wwTopoDesignON) return;
  
    //Map-based visualization transformation
    SoTransform *mapToWorldTrans = new SoTransform;
    mapSpaceProjector.getTransform(mapToWorldTrans);
    wwdvScene->addChild(mapToWorldTrans);
    
    //Picked Point: converted to a map location constrained to DeltaV line
    SbVec3f mapPoint = mapSpaceProjector.worldToMap(
                        horzMenu,vertMenu,horzScale,vertScale,
                        pickedPoint);
  
    //Augmented connection is initiated externally, so the data is just there...
    State &state0 = augHConn.manifoldState0;
    SbVec3f designStart(state0.x,state0.xd,0.0);
    
    //Render the constraint line (dotted, underneath)
    float constraintVerts[2][3] = {
        {designStart[0], portYbounds.getValue(0), 0.0},
        {designStart[0], portYbounds.getValue(1), 0.0}
    };  
    SoSeparator *lsSep = new SoSeparator;
    SoCoordinate3 *lsCoords = new SoCoordinate3;
    SoLineSet *lsLine = new SoLineSet;
    //Draw Style
    SoDrawStyle *lStyle = new SoDrawStyle;
    lStyle->lineWidth = portLineWidth.getValue();
    lStyle->linePattern.setValue( 10110 );
    SoMaterialBinding *lsMatBind = new SoMaterialBinding;
    //One color for each line segment
    lsMatBind->value = SoMaterialBinding::PER_PART;
    SoMaterial *lsMat = new SoMaterial;
    lsCoords->point.setValues(0,2,constraintVerts);
    lsLine->numVertices.set1Value(0,2);
    McColor cLineCol; designColors.getColor(cLineCol,2);
    lsMat->diffuseColor.set1Value(0,SbColor(cLineCol.r,cLineCol.g,cLineCol.b));
    //Assign design selection node
    topoDesignTargetNode = lsSep;
    //Add to display
    lsSep->addChild(lsMatBind);
    lsSep->addChild(lsMat);
    lsSep->addChild(lStyle);
    lsSep->addChild(lsCoords);
    lsSep->addChild(lsLine);
    //Insert to scene at fixed location (after event callback)
    wwdvScene->addChild(lsSep);
    
    //Render line from the starting point location to the pickedPoint (constrained by selection)
    SbVec3f *markerCoord = new SbVec3f[2];
    markerCoord[0] = SbVec3f(state0.x,state0.xd,0.0);  //Known state
    markerCoord[1] = SbVec3f(state0.x,mapPoint[1],0.0); //Constrained selection
    int markerIndex[2] = { 
        SoMarkerSet::CIRCLE_FILLED_9_9,
        SoMarkerSet::SATELLITE_FILLED_9_9
    };
    SbColor *markerColor = new SbColor[2];
    McColor iStateCol; designColors.getColor(iStateCol,0);
    markerColor[0].setValue(iStateCol.r,iStateCol.g,iStateCol.b);
    McColor ppCol; designColors.getColor(ppCol,1);
    //Only color working design with correct color:
    if (topoDesignSelected) {
        markerColor[1].setValue(ppCol.r,ppCol.g,ppCol.b);
    } else {
        markerColor[1].setValue(1.0,0.0,0.0);
    }
    SoCoordinate3 *markerCoord3 = new SoCoordinate3;
    markerCoord3->point.setValues(0,2,markerCoord);
    SoMarkerSet *markerSet = new SoMarkerSet;
    markerSet->markerIndex.setValues(0,2,markerIndex);
    markerSet->markerScale.setValue( portMarkerScale.getValue() );
    SoMaterialBinding *matBind = new SoMaterialBinding;
    matBind->value = SoMaterialBinding::PER_VERTEX;
    SoMaterial *markerMat = new SoMaterial;
    markerMat->diffuseColor.setValues(0,2,markerColor);
    wwdvScene->addChild(markerCoord3);
    wwdvScene->addChild(matBind);
    wwdvScene->addChild(markerMat);
    wwdvScene->addChild(markerSet);
      
    //Make sure initial state is within bounds selected by object
      // ->Skip for now...  assume ok...
    
    //Display the translation from an input point as vector 
    float dvVerts[2][3] = {
        {designStart[0], designStart[1], 0.0},
        {designStart[0], mapPoint[1], 0.0}
    };  
    SoSeparator *lsSep2 = new SoSeparator;
    SoCoordinate3 *lsCoords2 = new SoCoordinate3;
    SoLineSet *lsLine2 = new SoLineSet;    
    //Draw Style
    SoDrawStyle *lsStyle2 = new SoDrawStyle;
    lsStyle2->lineWidth = 2.0*portLineWidth.getValue();
    SoMaterialBinding *lsMatBind2 = new SoMaterialBinding;
    //One color for each line segment
    lsMatBind2->value = SoMaterialBinding::PER_PART;
    SoMaterial *lsMat2 = new SoMaterial;
    lsCoords2->point.setValues(0,2,dvVerts);
    lsLine2->numVertices.set1Value(0,2);
    McColor dvLineCol; designColors.getColor(dvLineCol,1);
    lsMat2->diffuseColor.set1Value(0,SbColor(dvLineCol.r,dvLineCol.g,dvLineCol.b));
    lsSep2->addChild(lsMatBind2);
    lsSep2->addChild(lsMat2);
    lsSep2->addChild(lsStyle2);
    lsSep2->addChild(lsCoords2);
    lsSep2->addChild(lsLine2);
    wwdvScene->addChild(lsSep2);
    
    //Display Orbit ID? -> Rendered in selected form
    
    //Add the eventCB into the scene to look for specified topology node
    //wwdvScene->insertChild(eventCB,0);
    // eventCB should be at root level!
  
    showGeom(wwdvScene);
}

///Output selected manifold arc as HxLineSet
void DisplayManifolds::exportArcToLineSet(HxLineSet *line)
{
  //First, determine how to construct states (up,down,both)
  int numSteps = (int) uStates.size();
  bool upstreamSim = true, downstreamSim = false;
  switch (simDirection.getValue() ) {
    case 1 :
      upstreamSim = false; downstreamSim = true;
      numSteps = (int) dStates.size();
      break;
    case 2 :
      upstreamSim = true; downstreamSim = true;
      numSteps = (int) uStates.size() + (int) dStates.size();
      break;
    default :
      upstreamSim = true; downstreamSim = false;
      numSteps = (int) uStates.size();
      break;
  }
  
  //Construct the data containers
  McVec3f *vertex = new McVec3f[numSteps];
  int *Ind = new int[numSteps];
  float *val0 = new float[numSteps];
  float *val1 = new float[numSteps];
  float *val2 = new float[numSteps];
  float *val3 = new float[numSteps];
  float *val4 = new float[numSteps];

  //Populate data based on simulation choices
  if (upstreamSim) {
    for (int i=0; i<(int)uStates.size(); i++) {
          //Store values
          vertex[i].setValue(uStates[i].x,uStates[i].y,uStates[i].z);
          Ind[i] = i;
          val0[i] = uStates[i].xd;
          val1[i] = uStates[i].yd;
          val2[i] = uStates[i].zd;
          val3[i] = uTimes[i];
          val4[i] = (float) sys->getJacobiConstant(uStates[i]); /*Jacobi Constant*/
    }
  }
  if (downstreamSim) { //Tack on to upstream if both are on:
    for (int i=0; i<(int)dStates.size(); i++) {
          //Store values
          vertex[i].setValue(dStates[i].x,dStates[i].y,dStates[i].z);
          Ind[i] = i;
          val0[i] = dStates[i].xd;
          val1[i] = dStates[i].yd;
          val2[i] = dStates[i].zd;
          val3[i] = dTimes[i];
          val4[i] = (float) sys->getJacobiConstant(dStates[i]); /*Jacobi Constant*/
    }
  }

  //Add points to lineset
  line->addPoints(vertex,numSteps);
  delete vertex;
  line->addLine(numSteps,Ind);
  delete Ind;
  //Add Data Values
  line->setNumDataValues(5);
  memcpy(line->getData(0),val0,numSteps*sizeof(float));
  memcpy(line->getData(1),val1,numSteps*sizeof(float));
  memcpy(line->getData(2),val2,numSteps*sizeof(float));
  memcpy(line->getData(3),val3,numSteps*sizeof(float));
  memcpy(line->getData(4),val4,numSteps*sizeof(float));
  delete[] val0; delete[] val1; delete[] val2; delete[] val3; delete[] val4;
}

///Output selected manifold arc as HxTrajectory
void DisplayManifolds::exportArcToTrajectory(HxTrajectory *traj)
{
  //First, determine how to construct states (up,down,both)
  //int numSteps = (int) uStates.size();
  bool upstreamSim = true, downstreamSim = false;
  switch (simDirection.getValue() ) {
    case 1 :
      upstreamSim = false; downstreamSim = true;
      //numSteps = (int) dStates.size();
      break;
    case 2 :
      upstreamSim = true; downstreamSim = true;
      //numSteps = (int) uStates.size() + (int) dStates.size();
      break;
    default :
      upstreamSim = true; downstreamSim = false;
      //numSteps = (int) uStates.size();
      break;
  }
  
  //Construct the data containers
  std::vector<State> totalStates;
  std::vector<double> totalTimes;
  double elapsedTime = 0.0;
  
  //Store in same time direction (stable - backward time, unstable - forward time)
  if (upstreamSim) {
    for (int i=0; i<(int)uStates.size(); i++) {
          //Store values
          totalStates.push_back( uStates[i] );
          //Times
          totalTimes.push_back( uTimes[i] );
          elapsedTime = uTimes[i]; //Negative for stable manifold
    }
  }
  if (downstreamSim) { //Tack on to upstream if both are on:
    for (int i=0; i<(int)dStates.size(); i++) {
          //Store values
          totalStates.push_back( dStates[i] );
          //Times
          totalTimes.push_back( dTimes[i] + elapsedTime );
    }
  }

  //Commit to trajectory object 
  traj->data.setStates(totalStates, totalTimes);
  
  //Update object
  traj->data.name = "Selected manifold arc";
  traj->updatePorts();
  traj->compute();
  
  //Note: there are likely some discrepancies within this arc! 
  //Perhaps we need another class to handle and fix those small errors?
  // ->mypackage/MultipleShooting can fix the problem! Now has HxTrajectory Input!
}

///Output selected manifold returns to a HxCluster object
void DisplayManifolds::exportReturnsToCluster(HxCluster *cluster) 
{
    //Check for input data
    ManifoldData *mData = (ManifoldData*) portData.source();
    if(!mData) {   return;   }
    
    
    //Check for Connected FixedPointData
    FixedPointData *fpData = (FixedPointData*) mData->fpxDataConnection.source();
    if(!fpData) {  return;   }
    
    
    //Need some orbit data
    
    //First, determine how to construct states (up,down,both)
    int numSteps = (int) uIterates.size();
    bool upstreamSim = true, downstreamSim = false;
    switch (simDirection.getValue() ) {
      case 1 :
        upstreamSim = false; downstreamSim = true;
        numSteps = (int) dIterates.size();
        break;
      case 2 :
        upstreamSim = true; downstreamSim = true;
        numSteps = (int) uIterates.size() + (int) dIterates.size();
        break;
      default :
        upstreamSim = true; downstreamSim = false;
        numSteps = (int) uIterates.size();
        break;
    }
    
    //Convert segment, line index, and parameter into manifold and state
    LineSegPair pickedManifold( selectedLineID, selectedPartID );
    LineSegPair hcPartnerManifold( hcPartnerLineID, hcPartnerPartID );
    //LineID,SegID -> ManifoldID
    IntPair manSegPair;
    IntPair hcManSegPair;
    if(selectedUnstable) {
      //Get (manifoldID,segmentID)
      manSegPair = uLinesetToManifold[pickedManifold];
    } else {
      manSegPair = sLinesetToManifold[pickedManifold];
    }
    //Heteroclinic piece
    if(hcOutput.getValue(0)==1) {
      if(selectedUnstable) {
        //Get (manifoldID,segmentID) for HC partner
        hcManSegPair = sLinesetToManifold[hcPartnerManifold];
      } else {
        hcManSegPair = uLinesetToManifold[hcPartnerManifold];
      }
    }
    
    //Get the Manifold and Segment 
    int mID = manSegPair.first; //int segID = manSegPair.second;
    const ManifoldData::ManifoldType& theManifold = mData->theManifoldData.mapManifolds[mID];
    //const ManifoldData::ManifoldSeg& theSeg = theManifold.segments[segID];
    int thePeriod = theManifold.getThePeriod();
    int orbitID = theManifold.fpdOrbitIdx;
    //int currentDepth = theManifold
    
  
    //Setup the HxCluster object for printing the returns
    double elapsedTime = 0.0;
    double elapsedReturns = 0;
    int numPts = numSteps;
    cluster->setNumDataColumns(6);
    cluster->setNumLabelColumns(1);
    cluster->resize(numPts);
    //Data Arrays
    //float *vMag = (float*) cluster->dataColumns[0].data;
    //cluster->dataColumns[0].name = "Velocity Magnitude";
    //float *ydotVal = (float*) cluster->dataColumns[1].data;
    //cluster->dataColumns[1].name = "Velocity ydot";
    float *tVal = (float*) cluster->dataColumns[0].data;
    cluster->dataColumns[0].name = "Time (NonDim)";
    float *crossVal = (float*) cluster->dataColumns[1].data;
    cluster->dataColumns[1].name = "Crossing Number";
    float *jcVal = (float*) cluster->dataColumns[2].data;
    cluster->dataColumns[2].name = "Jacobi";
    float *oID = (float*) cluster->dataColumns[3].data;
    cluster->dataColumns[3].name = "OrbitID";
    float *manID = (float*) cluster->dataColumns[4].data;
    cluster->dataColumns[4].name = "ManID";
    float *fwd = (float*) cluster->dataColumns[5].data;
    cluster->dataColumns[5].name = "Fwd";
    //float *dLevel = (float*) cluster->dataColumns[5].data;
    //cluster->dataColumns[5].name = "Depth";
    McString *crossString = (McString*) cluster->labelColumns[0].data;
    cluster->labelColumns[0].name = "Return Label";
    
    //Store in same time direction (stable - backward time, unstable - forward time)
    if (upstreamSim) {
      for (int n=0; n<(int)uIterates.size(); n++) {
           //Store values
           //Have to convert states from "map" space to "world" space (i.e., display coords for the map)
           SbVec3f mapPoint(uIterates[n][0],uIterates[n][1],0.0);
           SbVec3f wPoint = mapSpaceProjector.mapToWorld(horzMenu,vertMenu,horzScale,vertScale,mapPoint);
           cluster->points[n].setValue(wPoint[0],wPoint[1],wPoint[2]);
           
           //Date values
           //vMag[n] = (float) uIterates[n].vmag();
           //ydotVal[n] = (float) iterates[n].yd;
           crossVal[n] = (float) n*thePeriod;
           elapsedReturns = n*thePeriod;
           QString nString; nString.setNum(n);
           crossString[n] = nString.toLocal8Bit().data();
           jcVal[n] = (float) fpData->getJacobiConstant();
           oID[n] = (float) orbitID;
           manID[n] = (float) mID;
           //dLevel[n] = (float) what;
           fwd[n] = (float) selectedUnstable;
           //Times
           tVal[n] = (float) uItTimes[n];
           elapsedTime = uItTimes[n]; //Negative for stable manifold
      }
    }
    if (downstreamSim) { //Tack on to upstream if both are on:
      for (int n=0; n<(int)dIterates.size(); n++) {
           //Have to convert states from "map" space to "world" space (i.e., display coords for the map)
           SbVec3f mapPoint(dIterates[n][0],dIterates[n][1],0.0);
           SbVec3f wPoint = mapSpaceProjector.mapToWorld(horzMenu,vertMenu,horzScale,vertScale,mapPoint);
           cluster->points[n].setValue(wPoint[0],wPoint[1],wPoint[2]);
           
           //Date values
           crossVal[n] = (float) (n*thePeriod+elapsedReturns);
           QString nString; nString.setNum(n);
           crossString[n] = nString.toLocal8Bit().data();
           jcVal[n] = (float) fpData->getJacobiConstant();
           oID[n] = (float) orbitID;
           manID[n] = (float) mID;
           //dLevel[n] = (float) what;
           fwd[n] = (float) selectedUnstable;
           //Times
           tVal[n] = (float) (dItTimes[n] + elapsedTime);
      }
    }
    
    //Register
    theObjectPool->addObject(cluster);
}

/// Print (and sort/unique) the orbitID list to the list text port
void DisplayManifolds::printListToPort() {
    //Sort list
    orbitIDList.sort();
    //Make unique
    orbitIDList.unique();
    //Create the string 
    QString idListStr;
    std::list<int>::iterator lit;
    for(lit=orbitIDList.begin();lit!=orbitIDList.end();lit++) {
      QString eleStr;
      eleStr.setNum((*lit));
      idListStr += eleStr + ", ";
    }
    //Set the port
    portOrbitIDList.setSensitivity(true);
    portOrbitIDList.setValue(idListStr);
    portOrbitIDList.setSensitivity(false);
}

/// Compute the topology design DeltaV vector (after selection and only when Initial State input available)
bool DisplayManifolds::topologyDeltaV()
{
    //Input state : mark map location (x,xdot)
    HxInitialState* iState = (HxInitialState*) initStateConnection.source();
    
    //InitialState data
    State &theState = iState->data.ic;
    double inputC = sys->getJacobiConstant(theState);
    
    //Make sure that we have selected something to trigger maneuver computation
    //if (!topoDesignSelected) {
    //  //Set to zero 
    //  deltaV.setValue(0.0,0.0,0.0);
    //  return;
    //}
    
    //Manifold and Fixed point data (will be available)
    ManifoldData *mData = (ManifoldData*) portData.source();
    FixedPointData *fpData = (FixedPointData*) mData->fpxDataConnection.source();
    
    
    //CR3BP info
    double mup = sys->mup;
    double jcManifolds = fpData->getJacobiConstant();
    rhs_type   theEOMs(jcManifolds,mup);  //The right-hand side:  CR3BP EoMs
    section_type  theSection(theEOMs); //Hyperplane: y=0 by planar_section
    
    
    ///Convert segment, line index, and parameter into manifold and state
    LineSegPair pickedManifold( selectedLineID, selectedPartID );
    //LineSegPair hcPartnerManifold( hcPartnerLineID, hcPartnerPartID );
    //LineID,SegID -> ManifoldID
    IntPair manSegPair;
    //IntPair hcManSegPair;
    if(selectedUnstable) {
      //Get (manifoldID,segmentID)
      manSegPair = uLinesetToManifold[pickedManifold];
    } else {
      manSegPair = sLinesetToManifold[pickedManifold];
    }
    /*//Heteroclinic piece - Ignore
    if(hcOutput.getValue(0)==1) {
      if(selectedUnstable) {
        //Get (manifoldID,segmentID) for HC partner
        hcManSegPair = sLinesetToManifold[hcPartnerManifold];
      } else {
        hcManSegPair = uLinesetToManifold[hcPartnerManifold];
      }
    }*/
    
    //The picked point
    SbVec3f mapPoint; //(x,xdot,ydot)
    mapPoint = mapSpaceProjector.worldToMap(horzMenu,vertMenu,horzScale,vertScale,pickedPoint);
    nvis::vec2 xApprox(mapPoint[0],mapPoint[1]);
    //Get the Manifold and Segment 
    int mID = manSegPair.first; int segID = manSegPair.second;
    const ManifoldData::ManifoldType& theManifold = mData->theManifoldData.mapManifolds[mID];
    const ManifoldData::ManifoldSeg& theSeg = theManifold.segments[segID];
    //Compute the linear parameter for selected point (as a guess)
    selectedLinearParameter = theSeg.getLinearParam(xApprox);
    
    //Construct a dummy manifold segment to hold the deltaV-line information
    ManifoldData::ManifoldSeg dVLineSeg( 
         nvis::vec2(theState.x,(double)portYbounds.getValue(0)),
         nvis::vec2(theState.x,(double)portYbounds.getValue(1)),
         nvis::vec3(0.0) ); //Dummy data [not used]
    double dummyTau = 0.0;
    
    //Compute parameter based on intersection of deltaV line and the selected segment
    if ( !theSeg.intersect(dVLineSeg) ) {
      theMsg->printf("%s : Warning! DeltaV Line and selected manifold do not intersect!",getName());
    }
    bool found = theSeg.getIntersection(dVLineSeg,selectedLinearParameter,dummyTau);
    
    if(!found) return false;
    
    //The actual selected manifold point to double precision:
    nvis::vec2 xPP = theSeg.getPoint(selectedLinearParameter); //MapState (x,xdot)
    
    
    // Get the full state of the manifold with a reverse projection call
    State ppState;
    vec42 y = theSection.unproject(xPP);
    for(int i=0;i<6;i++) ppState[i] = y[i];

    //Choice on direction:  
    bool validPoint = true;
    double dxd = 0.0, dyd = 0.0, dC = 0.0, disc = 0.0;
    bool fwd = true;
    //Split settings based on design selection :
    //--------------------------------------------------------------------------------------
    //State to Ws (F) || State to Wu (B)
    if (portDesignOpt.getValue() == 0 || portDesignOpt.getValue() == 2) {
        //Set reverse design direction if going backward
        if(portDesignOpt.getValue() == 2) fwd = false;
        //Translation 
        dxd = xPP[1] - theState.xd;
        //Check sign!!!
        if (!fwd) dxd *= -1.0;
        //Jacobi difference
        dC = jcManifolds - inputC;
        //MapManeuver: InitialState to picked point
        disc = theState.yd*theState.yd - 2.0*theState.xd*dxd-dxd*dxd;
        //Only apply dC change when large enough
        if (std::fabs(dC)>= 1e-8) {
            //Free form definition
            disc = theState.yd*theState.yd - 2.0*theState.xd*dxd - dxd*dxd - dC;
        } 
        //Check zero-velocity condition, only proceed if valid 
        if (disc < 0.0) {
            //Set the selection point as red, indicating invalid!
            validPoint = false;
            //Keep values at zero.
        } else {
            //Compute the cheapest MAGNITUDE
            double dydot0 = -theState.yd + sqrt(disc);
            double dydot1 = -theState.yd - sqrt(disc);
            dyd = (std::fabs(dydot0)<std::fabs(dydot1))? dydot0 : dydot1;
        }
    }
    //--------------------------------------------------------------------------------------
    //Wu to State (F) || Ws to State (B)
    else if (portDesignOpt.getValue() == 1 || portDesignOpt.getValue() == 3) { 
        //Set reverse design direction if going backward
        if(portDesignOpt.getValue() == 3) fwd = false;
        //Translation
        dxd = theState.xd - xPP[1];
        //Check sign!!!
        if (!fwd) dxd *= -1.0;
        //Jacobi difference
        dC = inputC - jcManifolds;
        //MapManuever: PickedPoint to Input State
        disc = ppState.yd*ppState.yd - 2.0*ppState.xd*dxd - dxd*dxd;
        //Only apply dC change when large enough
        if (std::fabs(dC)>= 1e-8) {
            //Free form definition
            disc = ppState.yd*ppState.yd - 2.0*ppState.xd*dxd - dxd*dxd - dC;
        } 
        //Check zero-velocity condition, only proceed if valid 
        if (disc < 0.0) {
            //Set the selection point as red, indicating invalid!
            validPoint = false;
            //Keep values at zero.
        } else {
            //Compute the cheapest MAGNITUDE
            double dydot0 = -ppState.yd + sqrt(disc);
            double dydot1 = -ppState.yd - sqrt(disc);
            dyd = (std::fabs(dydot0)<std::fabs(dydot1))? dydot0 : dydot1;
        }
    }

    //These don't apply to InitialState input (hasInput = false)
    //--------------------------------------------------------------------------------------
       //Wu to Ws (F)
       //Ws to Wu (B)
    
    
    //Set the deltaV values in non-dim coords:
    if(validPoint) {
        topoDesignSelected = true;
        deltaV.setValue(dxd,dyd,0.0);
    } else {
        topoDesignSelected = false;
        //Set to zero 
        deltaV.setValue(0.0,0.0,0.0);
        //Prompt user 
        theMsg->printf("%s: Error!  Invalid maneuver based on selection (outside ZVC)!",__FILE__);
    }

    return topoDesignSelected;
}

/// Function to construct the heteroclinic/homoclinic data pairing (useful for design information)
bool DisplayManifolds::makeHHConnection()
{
  ManifoldData *mData = (ManifoldData*) portData.source();    
  //FixedPointData *fpData = (FixedPointData*) mData->fpxDataConnection.source();
  
  rhs_type theEOMs(C,mup);  //The right-hand side:  CR3BP EoMs
  section_type theSection(theEOMs); //Hyperplane: y=0 by planar_section
  vec42 y(0.0);
    
  
  //Gather selection - convert to ManifoldData indexes
  // ->Convert segment, line index, and parameter into manifold and state
  LineSegPair pickedManifold( selectedLineID, selectedPartID );
  LineSegPair hcPartnerManifold( hcPartnerLineID, hcPartnerPartID );
  //LineID,SegID -> ManifoldID
  IntPair manSegPair;
  IntPair hcManSegPair;
  if(selectedUnstable) {
    //Get (manifoldID,segmentID)
    manSegPair = uLinesetToManifold[pickedManifold];
    //Get (manifoldID,segmentID) for HC partner
    hcManSegPair = sLinesetToManifold[hcPartnerManifold];
  } else {
    manSegPair = sLinesetToManifold[pickedManifold];
    hcManSegPair = uLinesetToManifold[hcPartnerManifold];
  }
    
  //Convert the picked point into map space with the projector
  SbVec3f mapPoint; //(x,xdot,ydot)
  mapPoint = mapSpaceProjector.worldToMap(horzMenu,vertMenu,horzScale,vertScale,pickedPoint);
  nvis::vec2 xApprox(mapPoint[0],mapPoint[1]);
  //Get the Manifold and Segment 
  int mID = manSegPair.first; int segID = manSegPair.second;
  const ManifoldData::ManifoldType& theManifold = mData->theManifoldData.mapManifolds[mID];
  const ManifoldData::ManifoldSeg& theSeg = theManifold.segments[segID];
  //Compute the linear parameter for selected point
  selectedLinearParameter = theSeg.getLinearParam(xApprox);
  //The partner
  int hcManID = hcManSegPair.first; int hcSegID = hcManSegPair.second;
  const ManifoldData::ManifoldType& theHCMan = mData->theManifoldData.mapManifolds[hcManID];
  const ManifoldData::ManifoldSeg& theHCSeg = theHCMan.segments[hcSegID];
  //Compute the linear parameter for selected point
  double hcLP = theHCSeg.getLinearParam(xApprox);
  
  //Construct intersection with ManifoldSegment class members : Performs the intersection of two segments
  /*bool checkIntersection = theSeg.intersect(theHCSeg); //Should be true
  if(!checkIntersection) {
      //Warn user
      theMsg->printf("%s : Heteroclinic intersection is detected but found no intersection!",getName());
  }
  printf(" Heteroclinic Intersection info: \n");
  printf(" --------------------------------------------------------------------------\n");
  printf("  slID = %d spID = %d   | hclID = %d hcpID = %d\n",
         selectedLineID, selectedPartID, hcPartnerLineID, hcPartnerPartID);
  
  printf("  Seg0 : [%g,%g] -> [%g,%g]\n",theSeg[0][0],theSeg[0][1],theSeg[1][0],theSeg[1][1]);
  printf("  Seg1 : [%g,%g] -> [%g,%g]\n",theHCSeg[0][0],theHCSeg[0][1],theHCSeg[1][0],theHCSeg[1][1]);
  printf("  Intersection Found = %s\n", ((checkIntersection)?"Yes":"NO") );
  printf("  Guess Params:  tau0 = %g  tau1 = %g\n", selectedLinearParameter, hcLP);//From xApprox*/
  bool found = theSeg.getIntersection(theHCSeg, selectedLinearParameter, hcLP);
  //printf("  Acutal Test: %s  tau0 = %g tau1 = %g\n",((found)?"FOUND":"NONE"),selectedLinearParameter,hcLP);
  
  //It's possible that the "detected" intersection is not actually an intersection
  if(!found) return false;
  
  
  //Store data into design object [ALWAYS Wu to Ws]
  hhConn.fwd = true;
  if(selectedUnstable) {
    //Commons
    hhConn.tof = theManifold.getTimeOfFlight(segID,selectedLinearParameter);
    hhConn.tof += std::fabs( theHCMan.getTimeOfFlight(hcSegID,hcLP) );
    hhConn.isValid = found;
    //The unstable manifold arc
    hhConn.manID0 = mID;
    hhConn.segID0 = segID;
    hhConn.tau0 = selectedLinearParameter;
    hhConn.orbitID0 = theManifold.fpdOrbitIdx;
    hhConn.fpID0 = theManifold.fpdPointIdx;
    hhConn.man0LineIdx = selectedLineID;
    hhConn.man0PartIdx = selectedPartID;
    y = theSection.unproject(theSeg.getPoint(selectedLinearParameter));
    for(int i=0;i<6;i++) hhConn.manifoldState0[i] = y[i];
    for(int i=0;i<3;i++) {
        hhConn.deltaV[i] = y[i]; //Position of transfer point (On-section)
        hhConn.deltaV[i+3] = 0.0; //Free transfer 
    }
    //The stable manifold arc 
    hhConn.manID1 = hcManID;
    hhConn.segID1 = hcSegID;
    hhConn.tau1 = hcLP;
    hhConn.orbitID1 = theHCMan.fpdOrbitIdx;
    hhConn.fpID1 = theHCMan.fpdPointIdx;
    hhConn.man1LineIdx = hcPartnerLineID;
    hhConn.man1PartIdx = hcPartnerPartID;
    y = theSection.unproject(theHCSeg.getPoint(hcLP));
    for(int i=0;i<6;i++) hhConn.manifoldState1[i] = y[i];
    
  } else {
    //Commons
    hhConn.tof = std::fabs( theManifold.getTimeOfFlight(segID,selectedLinearParameter) );
    hhConn.tof += theHCMan.getTimeOfFlight(hcSegID,hcLP);
    hhConn.isValid = found;
    //The stable manifold arc
    hhConn.manID1 = mID;
    hhConn.segID1 = segID;
    hhConn.tau1 = selectedLinearParameter;
    hhConn.orbitID1 = theManifold.fpdOrbitIdx;
    hhConn.fpID1 = theManifold.fpdPointIdx;
    hhConn.man1LineIdx = selectedLineID;
    hhConn.man1PartIdx = selectedPartID;
    y = theSection.unproject(theSeg.getPoint(selectedLinearParameter));
    for(int i=0;i<6;i++) hhConn.manifoldState1[i] = y[i];
    for(int i=0;i<3;i++) {
        hhConn.deltaV[i] = y[i]; //Position of transfer point (On-section)
        hhConn.deltaV[i+3] = 0.0; //Free transfer 
    }
    //The unstable manifold arc 
    hhConn.manID0 = hcManID;
    hhConn.segID0 = hcSegID;
    hhConn.tau0 = hcLP;
    hhConn.orbitID0 = theHCMan.fpdOrbitIdx;
    hhConn.fpID0 = theHCMan.fpdPointIdx;
    hhConn.man0LineIdx = hcPartnerLineID;
    hhConn.man0PartIdx = hcPartnerPartID;
    y = theSection.unproject(theHCSeg.getPoint(hcLP));
    for(int i=0;i<6;i++) hhConn.manifoldState0[i] = y[i];
  }
  
  return true;
}

/// Initiating a Manifold-to-Manifold Design by storing information with starting selection
void DisplayManifolds::startAugmentedConnection()
{
    //----------------------------------------------------------------------------------------
    //Gather information about the selection for WtoW design
    //----------------------------------------------------------------------------------------
    //Input IS available to start clicking
    ManifoldData *mData = (ManifoldData*) portData.source();
    //FixedPointData *fpData = (FixedPointData*) mData->fpxDataConnection.source();
    
    //Convert segment, line index, and parameter into manifold and state
    LineSegPair pickedManifold( selectedLineID, selectedPartID );
    //LineSegPair hcPartnerManifold( hcPartnerLineID, hcPartnerPartID );
    //LineID,SegID -> ManifoldID
    IntPair manSegPair;
    //IntPair hcManSegPair;
    if(selectedUnstable) {
      //Get (manifoldID,segmentID)
      manSegPair = uLinesetToManifold[pickedManifold];
    } else {
      manSegPair = sLinesetToManifold[pickedManifold];
    }
    //Heteroclinic piece - Ignore for WtoW design
    /*if(hcOutput.getValue(0)==1) {
      if(selectedUnstable) {
        //Get (manifoldID,segmentID) for HC partner
        hcManSegPair = sLinesetToManifold[hcPartnerManifold];
      } else {
        hcManSegPair = uLinesetToManifold[hcPartnerManifold];
      }
    }*/
    
    //Convert the picked point into map space with the projector
    SbVec3f mapPoint; //(x,xdot,ydot)
    mapPoint = mapSpaceProjector.worldToMap(horzMenu,vertMenu,horzScale,vertScale,pickedPoint);
    nvis::vec2 xApprox(mapPoint[0],mapPoint[1]);
    //Get the Manifold and Segment 
    int mID = manSegPair.first; int segID = manSegPair.second;
    const ManifoldData::ManifoldType& theManifold = mData->theManifoldData.mapManifolds[mID];
    const ManifoldData::ManifoldSeg& theSeg = theManifold.segments[segID];
    
    //Compute the linear parameter for selected point
    selectedLinearParameter = theSeg.getLinearParam(xApprox);
    
    //Utilizing a linear parameter is more accurate for manifold seeding 
    //than just a point (float) picked from manifold
    //    ManifoldID,segmentID,linParam -> map state (x,xdot) 
    nvis::vec2 x0 = theSeg.getPoint(selectedLinearParameter);
    vec42 y(0.0);
    rhs_type theEOMs(C,mup);  //The right-hand side:  CR3BP EoMs
    section_type theSection(theEOMs); //Hyperplane: y=0 by planar_section
    y = theSection.unproject(x0);
    
    //----------------------------------------------------------------------------------------
    //Assign to WWDesign object
    //----------------------------------------------------------------------------------------
    augHConn.reset(); //Make sure we have a new set of connection information.
    //augHConn.isValid = false; //Not a good connection until we know it!
    augHConn.fwd = selectedUnstable;
    augHConn.manID0 = mID;
    augHConn.segID0 = segID;
    augHConn.tau0 = selectedLinearParameter;
    augHConn.orbitID0 = theManifold.fpdOrbitIdx;
    augHConn.fpID0 = theManifold.fpdPointIdx;
    augHConn.man0LineIdx = selectedLineID;
    augHConn.man0PartIdx = selectedPartID;
    for(int i=0;i<6;i++) augHConn.manifoldState0[i] = y[i];
    //Debug:
    //theMsg->printf("Start AugHConn: fwd = %d manID0 = %d segID0 = %d tau0 = %g orbitID0 = %d fpID0 = %d",
    //               selectedUnstable, mID, segID, selectedLinearParameter, augHConn.orbitID0, augHConn.fpID0);
    //theMsg->printf("    ManState0 = [ %g %g %g  |  %g %g %g ]",
    //               augHConn.manifoldState0[0], augHConn.manifoldState0[1], augHConn.manifoldState0[2],
    //               augHConn.manifoldState0[3], augHConn.manifoldState0[4], augHConn.manifoldState0[5]);
    //
    //----------------------------------------------------------------------------------------
    //Set the ports and other members according to what was selected
    //----------------------------------------------------------------------------------------
    /*if(selectedUnstable) {
      portDesignOpt.setValue(0,4);
    } else {
      portDesignOpt.setValue(0,5);
    }*/
    deltaV.setValue(0.0,0.0,0.0);
    wwTopoDesignON = true; //Starting the Manifold to Manifold design
    
}

/// Compute the Manifold-to-Manifold Design based on a valid 'hover-over' a secondary manifold
bool DisplayManifolds::computeAugmentedConnection()
{
    //Manifold and Fixed point data (will be avaliable)
    ManifoldData *mData = (ManifoldData*) portData.source();
    FixedPointData *fpData = (FixedPointData*) mData->fpxDataConnection.source();
    
    
    //CR3BP info
    double mup = sys->mup;
    double jcManifolds = fpData->getJacobiConstant();
    rhs_type   theEOMs(jcManifolds,mup);  //The right-hand side:  CR3BP EoMs
    section_type  theSection(theEOMs); //Hyperplane: y=0 by planar_section
    
    //Other manifold state is pre-computed
    State &st0 = augHConn.manifoldState0;
    const ManifoldData::ManifoldType& theOther = mData->theManifoldData.mapManifolds[augHConn.manID0];
    double partialTof = theOther.getTimeOfFlight(augHConn.segID0,augHConn.tau0);
    
    ///Convert segment, line index, and parameter into manifold and state
    LineSegPair pickedManifold( selectedLineID, selectedPartID );
    //LineSegPair hcPartnerManifold( hcPartnerLineID, hcPartnerPartID );
    //LineID,SegID -> ManifoldID
    IntPair manSegPair;
    //IntPair hcManSegPair;
    if(selectedUnstable) {
      //Get (manifoldID,segmentID)
      manSegPair = uLinesetToManifold[pickedManifold];
    } else {
      manSegPair = sLinesetToManifold[pickedManifold];
    }
    /*//Heteroclinic piece - Ignore
    if(hcOutput.getValue(0)==1) {
      if(selectedUnstable) {
        //Get (manifoldID,segmentID) for HC partner
        hcManSegPair = sLinesetToManifold[hcPartnerManifold];
      } else {
        hcManSegPair = uLinesetToManifold[hcPartnerManifold];
      }
    }*/
    
    //Get the Manifold and Segment 
    int mID = manSegPair.first; int segID = manSegPair.second;
    const ManifoldData::ManifoldType& theManifold = mData->theManifoldData.mapManifolds[mID];
    const ManifoldData::ManifoldSeg& theSeg = theManifold.segments[segID];
    //Guess at Linear Param from approximate position
    SbVec3f mapPoint; //(x,xdot,ydot)
    mapPoint = mapSpaceProjector.worldToMap(horzMenu,vertMenu,horzScale,vertScale,pickedPoint);
    nvis::vec2 xApprox(mapPoint[0],mapPoint[1]);
    selectedLinearParameter = theSeg.getLinearParam(xApprox);
    
    //Construct a dummy manifold segment to hold the deltaV-line information
    ManifoldData::ManifoldSeg dVLineSeg( 
         nvis::vec2(st0.x,(double)portYbounds.getValue(0)),
         nvis::vec2(st0.x,(double)portYbounds.getValue(1)),
         nvis::vec3(0.0) ); //Dummy data [not used]
    double dummyTau = 0.0;
    
    //Compute parameter based on interesection of deltaV line and the selected segment
    /*if ( !theSeg.intersect(dVLineSeg) ) {
      theMsg->printf("%s : Warning! DeltaV Line and selected manifold do not intersect!",getName());
    }
    printf(" Augmented Connection Intersection info: \n");
    printf(" --------------------------------------------------------------------------\n");
    printf("  slID = %d spID = %d  \n", selectedLineID, selectedPartID);
    printf("  Seg0 : [%g,%g] -> [%g,%g]\n",theSeg[0][0],theSeg[0][1],theSeg[1][0],theSeg[1][1]);
    printf("  Seg1 : [%g,%g] -> [%g,%g]\n",dVLineSeg[0][0],dVLineSeg[0][1],dVLineSeg[1][0],dVLineSeg[1][1]);
    printf("  Intersection Found = %s\n", ((theSeg.intersect(dVLineSeg))?"Yes":"NO") );
    printf("  Guess Params:  tau = %g \n", selectedLinearParameter);//From xApprox*/
    bool found = theSeg.getIntersection(dVLineSeg,selectedLinearParameter,dummyTau);
    //printf("  Acutal Test: %s  tau = %g  dummyTau = %g\n",((found)?"FOUND":"NONE"),selectedLinearParameter,dummyTau);
    
    if(!found) return false;
    
    //The actual selected manifold point to double precision:
    nvis::vec2 xPP = theSeg.getPoint(selectedLinearParameter); //MapState (x,xdot)
    State ppState;
    vec42 y = theSection.unproject(xPP);
    for(int i=0;i<6;i++) ppState[i] = y[i];
    
    
    //Choose what to do based on design option:
    bool validPoint = true;
    double dxd = 0.0, dyd = 0.0, disc = 0.0;
    //----------------------------------------------------------------------------------------
    // Wu to Ws [Forward]
    if (portDesignOpt.getValue() == 4) {
      //Translation defines maneuver
      dxd = xPP[1] - st0.xd;
      //Note: At same Jacobi constant by design!
      //MapManeuver: Unstable Manifold State to Stable Manifold State
      disc = st0.yd*st0.yd - 2.0*st0.xd*dxd - dxd*dxd;
      //Check for zero-velocity condtion (shouldn't be an issue but just double-check)
      if (disc < 0.0) {
        //Set the selection point as red indicating invalid!
        validPoint = false;
        //Keep deltaV at zero
      } else {
          //Compute the cheapest MAGNITUDE
          double dydot0 = -st0.yd + sqrt(disc);
          double dydot1 = -st0.yd - sqrt(disc);
          dyd = (std::fabs(dydot0)<std::fabs(dydot1))? dydot0 : dydot1;
      }
    }
    //----------------------------------------------------------------------------------------
    // Ws to Wu [Backward]
    else {
      //Translation from Unstable(new selection) to Stable (st0)
      dxd = st0.xd - xPP[1]; //Forward (dt>0) deisgn
      //Note: At same Jacobi by design!
      disc = ppState.yd*ppState.yd - 2.0*ppState.xd*dxd - dxd*dxd;
      //Check for zero-velocity condtion (shouldn't be an issue but just double-check)
      if (disc < 0.0) {
        //Set the selection point as red indicating invalid!
        validPoint = false;
        //Keep deltaV at zero
      } else {
          //Compute the cheapest MAGNITUDE
          double dydot0 = -ppState.yd + sqrt(disc);
          double dydot1 = -ppState.yd - sqrt(disc);
          dyd = (std::fabs(dydot0)<std::fabs(dydot1))? dydot0 : dydot1;
      }
    }
    
    // Complete the Augmented Connection [Note: has to be valid point to be called]
    augHConn.manID1 = mID;
    augHConn.segID1 = segID;
    augHConn.tau1 = selectedLinearParameter;
    augHConn.manifoldState1 = ppState;
    augHConn.orbitID1 = theManifold.fpdOrbitIdx;
    augHConn.fpID1 = theManifold.fpdPointIdx;
    augHConn.man1LineIdx = selectedLineID;
    augHConn.man1PartIdx = selectedPartID;
    //Set common parameters
    if (validPoint) {
      augHConn.tof = std::fabs(partialTof) 
                   + std::fabs(theManifold.getTimeOfFlight(segID,selectedLinearParameter));
    }
    augHConn.deltaV = ppState; //Position
    augHConn.deltaV[3] = dxd;
    augHConn.deltaV[4] = dyd;
    augHConn.deltaV[5] = 0.0;
    deltaV.setValue(dxd,dyd,0.0);
    augHConn.isValid = validPoint; //Invalid if intersection not found
    
    return true;
}

/// Display the Manifold-to-Manifold topology design scene 
/// ->Only called when a valid WWDesign selection is made
/// ->Has option to clear
void DisplayManifolds::renderWWTopoDesignScene(bool remove)
{
    //Make this scene full of small objects as it is continuously updated
    tdScene->removeAllChildren();
    if(remove) return;
    
    //Render on top of stuff
    SoDepthOffset *tdOffset = new SoDepthOffset();
    tdScene->addChild(tdOffset);
    
    //Manifold and Fixed point data (will be avaliable)
    ManifoldData *mData = (ManifoldData*) portData.source();
    FixedPointData *fpData = (FixedPointData*) mData->fpxDataConnection.source();
    
 
    //CR3BP info
    double mup = sys->mup;
    double jcManifolds = fpData->getJacobiConstant();
    rhs_type   theEOMs(jcManifolds,mup);  //The right-hand side:  CR3BP EoMs
    section_type  theSection(theEOMs); //Hyperplane: y=0 by planar_section
    
    //Compute the DeltaV (Before rendering to check what's valid)
    bool validPoint = augHConn.isValid;
    if(!validPoint) return; //Nothing to render

    //Velocity parameters
    double dxd = deltaV[0], dvND = deltaV.length(); //NonDim
    
    //Get the Manifold and Segment for first object
    const ManifoldData::ManifoldType& theManifold0 = mData->theManifoldData.mapManifolds[augHConn.manID0];
    const ManifoldData::ManifoldSeg& theSeg0 = theManifold0.segments[augHConn.segID0];
    nvis::vec2 x0 = theSeg0.getPoint( augHConn.tau0 );
    SbVec3f mapPoint0(x0[0],x0[1],0.0);
    //Get the Manifold and Segment for second object
    const ManifoldData::ManifoldType& theManifold1 = mData->theManifoldData.mapManifolds[augHConn.manID1];
    const ManifoldData::ManifoldSeg& theSeg1 = theManifold1.segments[augHConn.segID1];
    nvis::vec2 x1 = theSeg1.getPoint( augHConn.tau1 );
    SbVec3f mapPoint1(x1[0],x1[1],0.0);
    

    //Transform 
    SoTransform *mapToWorldTrans = new SoTransform;
    mapSpaceProjector.getTransform(mapToWorldTrans);
    tdScene->addChild(mapToWorldTrans);

    //Display the starting and ending manifold points 
    SbVec3f *markerCoord = new SbVec3f[2];
    markerCoord[0] = SbVec3f(x0[0],x0[1],0.0);
    markerCoord[1] = SbVec3f(x1[0],x1[1],0.0);
    int markerIndex[2] = { 
        SoMarkerSet::CIRCLE_FILLED_9_9,
        SoMarkerSet::SATELLITE_FILLED_9_9
    };
    SbColor *markerColor = new SbColor[2];
    McColor iStateCol; designColors.getColor(iStateCol,0);
    markerColor[0].setValue(iStateCol.r,iStateCol.g,iStateCol.b);
    McColor ppCol; designColors.getColor(ppCol,1);
    //Working design (or no render)
    markerColor[1].setValue(ppCol.r,ppCol.g,ppCol.b);
    
    SoCoordinate3 *markerCoord3 = new SoCoordinate3;
    markerCoord3->point.setValues(0,2,markerCoord);
    SoMarkerSet *markerSet = new SoMarkerSet;
    markerSet->markerIndex.setValues(0,2,markerIndex);
    markerSet->markerScale.setValue( portMarkerScale.getValue() );
    SoMaterialBinding *matBind = new SoMaterialBinding;
    matBind->value = SoMaterialBinding::PER_VERTEX;
    SoMaterial *markerMat = new SoMaterial;
    markerMat->diffuseColor.setValues(0,2,markerColor);
    tdScene->addChild(markerCoord3);
    tdScene->addChild(matBind);
    tdScene->addChild(markerMat);
    tdScene->addChild(markerSet);
      
    //Make sure initial state is within bounds selected by object
      // ->Skip for now...  assume ok...
    
    //Display the translation (deltaV) vector 
    float dvVerts[2][3] = {
        {mapPoint0[0], mapPoint0[1], 0.0},
        {mapPoint1[0], mapPoint1[1], 0.0}
    };  
    SoSeparator *lsSep2 = new SoSeparator;
    SoCoordinate3 *lsCoords2 = new SoCoordinate3;
    SoLineSet *lsLine2 = new SoLineSet;    
    //Draw Style
    SoDrawStyle *lsStyle2 = new SoDrawStyle;
    lsStyle2->lineWidth = 2.0*portLineWidth.getValue();
    SoMaterialBinding *lsMatBind2 = new SoMaterialBinding;
    //One color for each line segment
    lsMatBind2->value = SoMaterialBinding::PER_PART;
    SoMaterial *lsMat2 = new SoMaterial;
    lsCoords2->point.setValues(0,2,dvVerts);
    lsLine2->numVertices.set1Value(0,2);
    McColor dvLineCol; designColors.getColor(dvLineCol,1);
    lsMat2->diffuseColor.set1Value(0,SbColor(dvLineCol.r,dvLineCol.g,dvLineCol.b));
    lsSep2->addChild(lsMatBind2);
    lsSep2->addChild(lsMat2);
    lsSep2->addChild(lsStyle2);
    lsSep2->addChild(lsCoords2);
    lsSep2->addChild(lsLine2);
    tdScene->addChild(lsSep2);

    //Render the deltaV text
    if (portDesignDisplayOpts.getValue(0) == 1) { 
      //Create a SoText2 to display the velocity magnitude
      SbVec3f dvTextMapCoord(mapPoint0[0],(mapPoint0[1]+mapPoint1[1])/2.0,0.0);
      SoMaterial *ppMat = new SoMaterial;
      McColor fontCol( portDesignFont.getFontColor() );
      ppMat->diffuseColor.setValue(fontCol.r,fontCol.g,fontCol.b);
      //Font
      SoFont *textFont = new SoFont;
      QString familyName = portDesignFont.getFontName();
      bool isBold = portDesignFont.isBoldFont();
      bool isItalic = portDesignFont.isItalicFont();
      if (isBold && isItalic) {
        familyName.append(" : Bold Italic");
      } else if (isBold) {
        familyName.append(" : Bold");
      } else if (isItalic) {
        familyName.append(" : Italic");
      }
      QByteArray fbarray = familyName.toLocal8Bit();
      textFont->name.setValue(fbarray.data());
      textFont->size.setValue(portDesignFont.getFontSize());
      //Add to sep to control all subsequent SoNodes
      SoSeparator *textSep = new SoSeparator;
      //textSep->addChild(textDOffset);
      textSep->addChild(ppMat);
      textSep->addChild(textFont);
      //Indicate a translation to define the location 
      SoTransform *textXform = new SoTransform;
      textXform->translation.setValue( dvTextMapCoord - SbVec3f(0,0,0) );
      //The actual Text 
      QString dvTextStr, dxdTextStr;
      double value = 0.0, value2 = 0.0;
      if (validPoint) {
          switch(portDeltaVUnits.getValue()) {
              case 1: //km/s 
                  //Convert the velocity magnitude to km/s 
                  value = sys->getDimVelocity(dvND);
                  dvTextStr.setNum(value,'g',5);
                  dvTextStr += " km/s";
                  value2 = sys->getDimVelocity(dxd);
                  dxdTextStr.setNum(value2,'g',5);
                  dxdTextStr += " km/s";
                  break;
              case 2: //m/s
                  //Convert the velocity magnitude to km/s 
                  value = sys->getDimVelocity(dvND);
                  dvTextStr.setNum(value*1.e3,'g',5);//m/s
                  dvTextStr += " m/s";
                  value2 = sys->getDimVelocity(dxd);
                  dxdTextStr.setNum(value2*1.e3,'g',5);
                  dxdTextStr += " m/s";
                  break;
              default : //ND 
                  dvTextStr.setNum(dvND,'g',5);
                  dvTextStr += " ND";
                  dxdTextStr.setNum(dxd,'g',5);
                  dxdTextStr += " ND";
                  break;
          }
      } else {
          dvTextStr = "N/A";
      }
      //Move off the line a bit and add vertical translation:
      QString fullText = " " + QString(QChar(0x0394)) + "V= " + dvTextStr;
      if (validPoint && portDesignDisplayOpts.getValue(2)==1) { //And another option?
        fullText += "  [" + QString(QChar(0x0394)) + "xd= " + dxdTextStr + "]";
      }
      SoText2 *dvText = new SoText2;
      SbString soStr; soStr.fromUtf16(fullText.utf16());
      dvText->string.setValue(soStr);
      //Add to scene 
      textSep->addChild(textXform);
      textSep->addChild(dvText);
      tdScene->addChild(textSep);
    }
    
    //Display the OrbitID of the selection
    if (portDesignDisplayOpts.getValue(1)==1) {
      //Create a SoText2 to display orbit ID of the selected manifold (second only!)
      
      //Get color from colormap
      SoMaterial *oidMat = new SoMaterial;
      SoMaterial *oidMat2 = new SoMaterial;
      if(augHConn.fwd) {
          oidMat->diffuseColor.setValue(unstableColormap.getColor(augHConn.orbitID0));
          oidMat2->diffuseColor.setValue(stableColormap.getColor(augHConn.orbitID1));
      } else {
          oidMat->diffuseColor.setValue(stableColormap.getColor(augHConn.orbitID0));
          oidMat2->diffuseColor.setValue(unstableColormap.getColor(augHConn.orbitID1));
      }
      //Font
      SoFont *textFont = new SoFont;
      QString familyName = portDesignFont.getFontName();
      bool isBold = portDesignFont.isBoldFont();
      bool isItalic = portDesignFont.isItalicFont();
      if (isBold && isItalic) {
        familyName.append(" : Bold Italic");
      } else if (isBold) {
        familyName.append(" : Bold");
      } else if (isItalic) {
        familyName.append(" : Italic");
      }
      QByteArray fbarray = familyName.toLocal8Bit();
      textFont->name.setValue(fbarray.data());
      textFont->size.setValue(portDesignFont.getFontSize());
      
      //Add to sep to control all subsequent SoNodes
      SoSeparator *textSep = new SoSeparator;
      //textSep->addChild(textDOffset);
      textSep->addChild(oidMat);
      textSep->addChild(textFont);
      //Indicate a translation to define the location 
      SoTransform *textXform = new SoTransform;
      textXform->translation.setValue( mapPoint0 - SbVec3f(0,0,0) );
      //The actual Text 
      QString oidTextStr;
      oidTextStr.setNum( augHConn.orbitID0 );
      //Move off the line a bit and add vertical translation:
      QString fullText = " ID = " + oidTextStr;
      SoText2 *oidText = new SoText2;
      oidText->string.setValue(fullText.toLocal8Bit().data());
      //The second point
      SoTransform *textXform2 = new SoTransform;
      textXform2->translation.setValue( mapPoint1 - mapPoint0 );
      QString oidTextStr2;
      oidTextStr2.setNum( augHConn.orbitID1 );
      QString fullText2 = " ID = " + oidTextStr2;
      SoText2 *oidText2 = new SoText2;
      oidText2->string.setValue(fullText2.toLocal8Bit().data());
      //Add to scene 
      textSep->addChild(textXform);
      textSep->addChild(oidText);
      textSep->addChild(oidMat2);
      textSep->addChild(textXform2);
      textSep->addChild(oidText2);
      tdScene->addChild(textSep);
    }
    
    
    //Show the scene 
    showGeom(tdScene);
}
