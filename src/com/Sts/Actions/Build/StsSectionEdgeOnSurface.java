
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Build;

import com.Sts.Actions.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;
import com.Sts.Interfaces.*;

import javax.swing.*;
import java.awt.event.*;

public class StsSectionEdgeOnSurface extends StsAction
{
    private StsModelSurface surface;
    private StsZone zone;
    private boolean visible = true;
    private StsSectionEdge edge;
    private StsSection section;
    private StsLine[] pickedLines = new StsLine[2];
    private int nLinesPicked = 0;
    private byte type;
    private String surfaceName = null;
	private boolean isNewSection = false;
	private StsObjectList availableSections;

	static private byte buttonType = StsParameters.FAULT;
	static private int nPointsPerVertex = 1;

   	/** create a new section and attach a mouse listener to it */
 	public StsSectionEdgeOnSurface(StsActionManager actionManager)
    {
        super(actionManager);
    }

	public void initializeRepeatAction(StsAction lastAction)
	{
		surface = ((StsSectionEdgeOnSurface)lastAction).getSurface();
	}

	public StsModelSurface getSurface() { return surface; }

	private boolean actionOK()
    {
    	try
        {
            StsModelSurfaceClass surfaceClass = (StsModelSurfaceClass)model.getCreateStsClass(StsModelSurface.class);
            if (surfaceClass.getSize() <= 0)
            {
                logMessage("No surfaces have been loaded: terminating action.");
                return false;
            }
	        return true;
        }
        catch(Exception e)
		{
			StsException.outputException("StsSectionEdgeOnSurface.actionOK() failed.",
				e, StsException.WARNING);
			return false;
		}
    }

	public boolean start()
    {
		if( !actionOK() )  return false;

        try
        {
            statusArea.setTitle("Build Fault Section: ");

            // select a horizon
			if(surface == null)
            {
                StsModelSurfaceClass surfaceClass = (StsModelSurfaceClass)model.getCreateStsClass(StsModelSurface.class);
                surface = (StsModelSurface)surfaceClass.selectSurface(true, StsModelSurface.class); // try visible 1st
            }
            if (surface == null)
            {
                statusArea.textOnly();
                logMessage("No horizon selected.");
                return false;
            }
            visible = surface.getIsVisible();
            surface.setIsVisible(true);

			surfaceName = surface.getName();

            if(!surface.checkIsLoaded())
            {
                statusArea.textOnly();
                logMessage("Unable to load " + surfaceName + ".");
                return false;
            }
            visible = surface.getIsVisible();
            surface.setIsVisible(true);


			logMessage("Select a fault line, digitize fault cut on " + surfaceName
                    + ", then finish by selecting another fault line.");

		    addAbortButton();
			addEndAllButton();
		    addEdgeTypeButton();

            /** Allocate a section; this may be temporary if we find
             * that we are adding an edge to an existing section */

//            model.viewProperties.save();
//            model.viewProperties.set("Sections", true);
//            model.viewProperties.set("Patches", false);

//            model.viewProperties.set("EditingSections", true);

            for(int n = 0; n < 2; n++) pickedLines[n] = null;

//            win3d.glPanel3d.set3dOverlay(true);
            model.win3d.win3dDisplay();
            return true;
        }
        catch(Exception e)
		{
			StsException.outputException("StsSectionEdgeOnSurface.start() failed.",
				e, StsException.WARNING);
			return false;
		}
    }

	public void addEdgeTypeButton()
    {
		String buttonText = StsParameters.typeToString(buttonType);
        statusArea.addActionButtonListener(new ButtonActionListener(), buttonText, true);
	}

    protected class ButtonActionListener implements ActionListener
    {
        public void actionPerformed(ActionEvent e)
        {
			String actionCommand = e.getActionCommand();
			if(actionCommand.equals("Fault"))
			{
				JButton actionButton = (JButton)e.getSource();
				actionButton.setActionCommand("Auxiliary");
				actionButton.setText("Auxiliary");
				buttonType = StsParameters.AUXILIARY;
			}
			else if(actionCommand.equals("Auxiliary"))
			{
				JButton actionButton = (JButton)e.getSource();
				actionButton.setActionCommand("Fault");
				actionButton.setText("Fault");
				buttonType = StsParameters.FAULT;
			}
        }
    }

    /** mouse action for 3d window */
   	public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel)
	{
        StsLine line;
        StsLine firstLinePicked = null;
		try
		{
            StsGLPanel3d glPanel3d = (StsGLPanel3d)glPanel;
            int leftButtonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);
			if(leftButtonState == StsMouse.PRESSED || leftButtonState == StsMouse.DRAGGED)
			{
				if(nLinesPicked == 1)
				{
					StsMousePoint mousePoint = mouse.getMousePoint();
					surface.getSurfacePosition(mouse, true, (StsGLPanel3d)glPanel);
				}
			}
			else if(StsMouse.isButtonStateReleasedOrClicked(leftButtonState))  // add a point
			{
				if(section == null) type = buttonType;

				/** first check if line has been picked; if so, add it to lineEnties list.
				  * If no line pick and one line has been picked,  compute intersection on surface
				  */

                StsObject[] visibleLines = getVisibleLines();
				line = (StsLine)StsJOGLPick.pickClass3d(glPanel3d, visibleLines, StsJOGLPick.PICKSIZE_MEDIUM, StsJOGLPick.PICK_CLOSEST);
				boolean lineVisible = line != null && line.getIsVisible();

				if (lineVisible && line != pickedLines[0])
				{
					while(line.hasSectionEdgeVertex(surface))
					{
					    logMessage("Fault/pseudo already has selected surface: select another.");
                        StsModelSurfaceClass surfaceClass = (StsModelSurfaceClass)model.getCreateStsClass(StsModelSurface.class);
						surface = (StsModelSurface)surfaceClass.selectSurface(false, StsModelSurface.class);
						if(surface == null) return true;
					}

					// maxConnections: number of sections which can be connected to this line
					// nConnectedSections: number of sections currently connected (connected to line or line is on section; max of 2)
					// nUnconnectedSections: maxConnections - nConnectedSections
					// nAvailableSections: number of allSections which do not have edge for this surface
					// nPossibleSections: nAvailableSections + nUnconnectedSections

					int maxConnections = line.getMaxConnections();
					StsSection[] allSections = line.getAllSections();
					int nConnectedSections = allSections.length;

					int nUnconnectedSections = maxConnections - nConnectedSections;
					StsSection[] availableSections = getAvailableSections(allSections);
					int nAvailableSections = availableSections.length;
					int nPossibleSections = nAvailableSections + nUnconnectedSections;

					if(nPossibleSections <= 0)
					{
						new StsMessage(model.win3d, StsMessage.WARNING, "No sections available to connect this fault/pseudo.");
						return true;
					}

					// if first line pick: highlight other lines on available sections
					if(nLinesPicked == 0)
					{
						if(nAvailableSections > 0)
						{
							for(int n = 0; n < nAvailableSections; n++)
							{
								StsSection availableSection = availableSections[n];
								StsLine otherLine = availableSection.getOtherEndLine(line);
								StsLine.highlightedList.add(otherLine);
							}
							if(nPossibleSections > nAvailableSections)
								logMessage("Picked " + line.getLabel() + ". Pick across " + surfaceName +
									". End on highlighted fault/pseudo or some other.");
							else if(nAvailableSections == 1)
							{
								logMessage("Picked " + line.getLabel() + ". Pick across " + surfaceName +
									". End on highlighted fault/pseudo.");
								section = availableSections[0];
							}
							else
								logMessage("Picked " + line.getLabel() + ". Pick across " + surfaceName +
									". End on one of the other two highlighted faults/pseudos.");
						}
						else // nAvailableSections == 0 && nPossibleSections > 0
							logMessage("Picked " + line.getLabel() + ". Pick across " + surfaceName +
								". End on another fault/pseudo.");

						StsLine.highlightedList.add(line);
					}
					pickedLines[nLinesPicked] = line;

//
//						line.addConnectedSection(section);

					// Need to improve this by computing nearest point on segment to
					// the mouse point, working in screen XY coordinates

					// find intersection of line line with grid, then snap to line line
					StsGridPoint gridPoint = line.computeGridIntersect(surface);
					StsPoint linePoint = line.getXYZPointAtZorT(gridPoint.getZorT(), true); // snap to line line

					nLinesPicked++;

					if(section == null)
					{
						isNewSection = true;
						// section = new StsSection(type);
					}
					StsSurfaceVertex vertex = new StsSurfaceVertex(linePoint, line, surface, null, true);
//					float dZMax = 2.0f * model.getProject().getZInc();
//                    line.insertSurfaceVertex(vertex);
//					line.insertLineVertex(vertex, false, 2*dZMax, true);

					if(nLinesPicked == 1)
					{
                        firstLinePicked = line;
                        // section.addEndLine(line);
					    edge = new StsSectionEdge(type, null, surface, nPointsPerVertex);
						model.addDisplayableInstance(edge);
						edge.addPrevVertex(vertex);
						edge.setDrawPoints(true);
					}
					else
					{
                        section = new StsSection(type, firstLinePicked, line);
                        // section.addEndLine(line);
						edge.addNextVertex(vertex);
						section.addSectionEdge(edge);

						if(isNewSection)
						{
						    StsSection duplicateSection = section.duplicateSection(availableSections);
						    if(duplicateSection != null)
						    {
							    duplicateSection.addSectionEdge(edge);
								section.delete();
								section = duplicateSection;
						    }
						}
						if(section.hasAllModelSurfaceEdges())
                        {
                            section.completeSection();
			                // section.constructSection();
						    glPanel.actionManager.endCurrentAction();
                        }

						edge.setDrawPoints(false);
					}
				}
				/** If line picked and it is same as firstLine, check screen Z values
				 *  to see if surface pick is behind line pick; if so, tell the user
				 *  that line has already been picked and skip out; if surface pick
				 *  is in front, than use it
				 */
				else if (line != null && line == pickedLines[0] && lineVisible)
				{
					StsPickItem pickItem = StsJOGLPick.pickItems[0];
					float lineScreenZ = pickItem.zMin;
					StsMousePoint mousePoint = mouse.getMousePoint();
					StsGridPoint gridPoint = surface.getSurfacePosition(mouse, true, (StsGLPanel3d)glPanel);
					double[] surfacePickCoor = glPanel3d.getScreenCoordinates(gridPoint.getPoint());
					if(lineScreenZ < (float)surfacePickCoor[2])
						new StsMessage(model.win3d, StsMessage.WARNING, "This line already picked.");
					else
					{
						StsGridSectionPoint edgePoint = new StsGridSectionPoint(new StsPoint(gridPoint.getPoint()), (StsSurfaceGridable)null, true);
						edge.addEdgePoint(edgePoint);
					}
				}
				 /** if no lines picked on this mouseAction, but one line has been picked,
				 * compute pick intersection with surface and add edgeVertex
				 */
				else if (nLinesPicked == 1)
				{
					StsMousePoint mousePoint = mouse.getMousePoint();
					StsGridPoint gridPoint = surface.getSurfacePosition(mouse, true, (StsGLPanel3d)glPanel);
					StsGridSectionPoint edgePoint = new StsGridSectionPoint(new StsPoint(gridPoint.getPoint()), (StsSurfaceGridable)null, true);
					edge.addEdgePoint(edgePoint);
				}
				else
					StsMessageFiles.logMessage("Must select a fault line before digitizing the surface cut.");
			}
			return true;
		}
		catch(Exception e)
		{
			StsException.outputException("StsSectionEdgeOnSurface.performMouseAction() failed.",
				e, StsException.WARNING);
			return false;
		}
		finally
		{

			if(edge != null)
			{
				edge.setPointsFromEdgePoints();
			    model.win3d.win3dDisplay();
			}
		}
	}

    private StsObject[] getVisibleLines()
    {
        StsObject[] visibleLines = model.getVisibleObjectList(StsLine.class);
        StsObject[] visibleFaultLines = model.getVisibleObjectList(StsFaultLine.class);
        return (StsObject[])StsMath.arrayAddArray(visibleLines, visibleFaultLines);
    }

	private StsSection[] getAvailableSections(StsSection[] connectedSections)
	{
        StsSection[] availableSections = new StsSection[0];
		if(connectedSections == null || connectedSections.length == 0) return availableSections;
		int nConnectedSections = connectedSections.length;
		for(int n = 0; n < nConnectedSections; n++)
		{
			if(!connectedSections[n].hasSurface(surface))
                availableSections = (StsSection[])StsMath.arrayAddElement(availableSections, connectedSections[n]);
		}
		return availableSections;
	}

  	public boolean end()
    {
        boolean ok = true;

		if(nLinesPicked == 0) return true; // nothing built; everthing ok
        surface.setIsVisible(visible);
		model.removeDisplayableInstance(edge);
//        model.viewProperties.restore();
		StsLine.clearHighlightedLines();
        statusArea.textOnly();

		if(nLinesPicked != 2)
		{
			if(isNewSection) section.delete();
			else edge.delete();
			ok = false;
            logMessage("Fault section not built.");
		}
		else
		{
            if(!section.initialized)
            {
                section.completeSection();
                // section.constructSection();
            }
            logMessage("Fault section:  " + section.getLabel() + " built successfully.");
			ok = true;
		}
//        win3d.glPanel3d.set3dOverlay(false);
        model.win3d.win3dDisplay();
        return ok;
	}
}


