
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Build;

import com.Sts.Actions.*;
import com.Sts.DBTypes.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.event.*;

public class StsEdgeOnSurface extends StsAction
{
    private StsXYSurfaceGridable surface; // current surface
	private StsXYSurfaceGridable[] surfaces; // array of top-down surfaces
	private int nSurfaces; // number of surfaces
	private boolean[] surfacePicked; // indicates surface already picked
	private int nEdges = 0;  // number of edges constructed
	private StsSurfaceVertex[][] edgeVertices; // [nEdges][nVertices] array
	private int surfaceIndex = 0; // curent index of edge
	private int nVerticesPicked = 0; // number of vertices picked for this edge
	private int firstVertexIndex = -1;  // vertex index in edgeVertices of first pick
	private int nVerticesMax = 0; // max nVertices determined from first edge picked
    private boolean visible = true;
    private StsSectionEdge edge;
    private StsObjectList sectionEdges = new StsObjectList(2, 2);
    private StsSection section;
	private int nSectionsBuilt = 0;
    private String surfaceName = null;
    private boolean ok = true;
	private int edgeDirection; // NONE, PLUS, or MINUS
//    private boolean addDepth = false;

    private byte type = StsParameters.FAULT;  // FAULT, or AUXILIARY
	static private int nPointsPerVertex = 1;

   	/** create a new edge on surface and attach a mouse listener to it */
 	public StsEdgeOnSurface(StsActionManager actionManager)
    {
        super(actionManager);
		setIsRepeatable(true);
		addUndoButton();
//		addEndAllButton();
    }

	public void initializeRepeatAction(StsAction lastAction)
	{
		surface = ((StsEdgeOnSurface)lastAction).getSurface();
	}

	public StsXYSurfaceGridable getSurface() { return surface; }

	public boolean start()
    {
		if( !actionOK() )  return false;

        try
        {
            statusArea.setTitle("Pick fault edges on surfaces: ");

            surfaces = (StsXYSurfaceGridable[])model.getCastObjectList(StsModelSurface.class);
            if (surfaces == null)
            {
                statusArea.textOnly();
                logMessage("No model surfaces available.");
                return false;
            }

//            model.viewProperties.save();

			initializeSurfaceDisplayFlags();
			initializeSectionPick();

			logMessage("Select a fault edge on surface or 'End Edge' button to skip.");

			addEndAllButton();
		    addActionButton();

            StsProject project = model.getProject();
            byte zDomainSupported = project.getZDomainSupported();
//            addDepth = supportedZDomain == StsProject.TD_TIME_DEPTH || supportedZDomain == StsProject.TD_APPROX_DEPTH_AND_DEPTH;

            model.win3d.win3dDisplay();
            return true;
        }
        catch(Exception e)
		{
			StsException.outputException("StsEdgeOnSurfaces.start() failed.",
				e, StsException.WARNING);
			return false;
		}
    }

	private void initializeSurfaceDisplayFlags()
	{
		nSurfaces = surfaces.length;
		for(int n = 0; n < nSurfaces; n++)
		{
			StsXYSurfaceGridable surface = surfaces[n];
            surface.toggleSurfacePickingOn();
            surface.setIsVisible(false);
        }
	}

	private void initializeSectionPick()
	{
		nSurfaces = surfaces.length;
		surfacePicked = new boolean[nSurfaces]; // assume initialized to false

		edgeVertices = new StsSurfaceVertex[0][];

		section = null;
	    nEdges = 0;
	    firstVertexIndex = -1;
	    nVerticesMax = 0;
		surfaceIndex = 0;

		setSurface();
	}

	private void setSurface()
	{
        if(surface != null)
            surface.toggleSurfacePickingOff();
		surface = surfaces[surfaceIndex];
        if(!surface.toggleSurfacePickingOn())
        {
            statusArea.textOnly();
            logMessage("Unable to load " + surfaceName + ".");
		    actionManager.endCurrentAction();
            return;
        }
        surface.setIsVisible(true);
        surfaceName = surface.getName();
	}

	public void addNextButton()
    {
        statusArea.addActionButtonListener(new ActionButtonListener(), "Next", true);
	}

	private boolean actionOK()
    {
    	try
        {
            StsObject[] surfaces = model.getObjectList(StsModelSurface.class);
            if(surfaces.length == 0)
            {
                logMessage("No surfaces have been loaded. Terminating action.");
                return false;
            }
	        return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsEdgeOnSurface.actionOK() failed.",
                e, StsException.WARNING);
            return false;
        }
    }

	public void addActionButton()
    {
		String buttonText = StsParameters.typeToString(type);
        statusArea.addActionButtonListener(new ActionButtonListener(), buttonText, true);
	}

    protected class ActionButtonListener implements ActionListener
    {
        public void actionPerformed(ActionEvent e)
        {
			String actionCommand = e.getActionCommand();
			if(actionCommand.equals("Fault"))
			{
				JButton actionButton = (JButton)e.getSource();
				setActionButtonCommand("Auxiliary");
				type = StsParameters.AUXILIARY;
			}
			else if(actionCommand.equals("Auxiliary"))
			{
				JButton actionButton = (JButton)e.getSource();
				setActionButtonCommand("Fault");
				type = StsParameters.FAULT;
			}
			else if(actionCommand.equals("End Edge"))
			{
				if(!completeEdge()) return;
				model.win3d.win3dDisplay();
//				if(section != null) setActionButtonCommand("End Section");
			}
			else if(actionCommand.equals("End Section"))
			{
				completeSection();
				model.win3d.win3dDisplay();
			}
        }
    }

	public void addUndoButton()
    {
        statusArea.addUndoButtonListener(new UndoButtonListener(), "Undo", true);
	}

    protected class UndoButtonListener implements ActionListener
    {
        public void actionPerformed(ActionEvent e)
        {
			if(edge != null) deleteEdge();
			else if(section != null) deleteSection();
        }
    }

	private void deleteEdge()
	{
		if(edge != null)
        {
            edge.delete();
            model.removeDisplayableInstance(edge);
        }
		if(edgeVertices != null) edgeVertices = (StsSurfaceVertex[][])StsMath.arrayDeleteLastElement(edgeVertices);
		nEdges = edgeVertices.length;
		if(nEdges == 0) initializeSectionPick();
		clearEdge();
        sectionEdges.delete(edge);
        model.win3d.win3dDisplay();
	}

    // mouse action for 3d window
   	public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel)
	{
		StsFaultLine pickedLine = null;
		StsSectionEdge pickedEdge = null;
		StsGridSectionPoint gridSectionPoint;
		float screenZLine = StsParameters.largeFloat;
		float screenZEdge = StsParameters.largeFloat;
		boolean edgePicked = false;
		StsPickItem pickItem;
		StsPoint pickedPoint = null;

		try
		{
            StsGLPanel3d glPanel3d = (StsGLPanel3d)glPanel;
            int leftButtonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);

			if(leftButtonState == StsMouse.PRESSED || leftButtonState == StsMouse.DRAGGED)
			{
				surface.getSurfacePosition(mouse, true, (StsGLPanel3d)glPanel);
				return true;
			}

			if(leftButtonState != StsMouse.RELEASED)  return true;

			// mouseButton RELEASED: process selection

			// haven't start edge yet: set action from Fault/Aux selection to "End Edge"
			if(edge == null) setActionButtonCommand("End Edge");

            // get intersection point on surface (if it exists)
            StsGridPoint gridPoint = surface.getSurfacePosition(mouse, true, (StsGLPanel3d)glPanel);
//            if(gridPoint == null) return true;

			// see if a line is picked
			pickedLine = (StsFaultLine)StsJOGLPick.pickClass3d(glPanel3d, (StsObject[])getConnectableLines(), StsJOGLPick.PICKSIZE_MEDIUM, StsJOGLPick.PICK_CLOSEST);
			if(pickedLine != null)
			{
				if(pickedLine.isFullyConnected())
				{
					StsMessage.printMessage("Sorry this fault is fully connected.\n" +
						"Select a point on an adjoining section or another fault.");
					pickedLine = null;
				}
				if(pickedLine != null)
				{
					pickItem = StsJOGLPick.pickItems[0];
					int nSegment = pickItem.names[1];
					pickedPoint = pickedLine.getPointOnLineNearestMouse(nSegment, mouse, glPanel3d);

                    if(gridPoint != null)
                    {
                        float[] xyz = pickedPoint.getXYZorT();
					    screenZLine = glPanel3d.getScreenZ(xyz);
//			            float screenZSurface = glPanel3d.getScreenZ(gridPoint.getXYZorT());
//				        if(screenZLine >= screenZSurface) pickedLine = null;
                    }
				}
			}

			pickedEdge = (StsSectionEdge)StsJOGLPick.pickClass3d( glPanel3d, getPickableEdges(), StsJOGLPick.PICKSIZE_MEDIUM, StsJOGLPick.PICK_CLOSEST);
			if(pickedEdge != null)
			{
				pickItem = StsJOGLPick.pickItems[0];
				int nSegment = pickItem.names[1];
				pickedPoint = pickedEdge.getPointOnLineNearestMouse(nSegment, mouse, glPanel3d);
//			    screenZEdge = win3d.glPanel3d.getScreenZ(pickedPoint.v);
//				if(screenZEdge >= screenZSurface) pickedEdge = null;
			}

			// if we have both a line and edge pick, delete the one furthest away
			if(pickedLine != null && pickedEdge != null)
			{
				if(screenZLine < screenZEdge)
					pickedEdge = null;
				else
					pickedLine = null;
			}

			if(pickedLine != null)
			{
				StsSurfaceVertex lineSurfaceVertex = pickedLine.getSectionEdgeVertex(surface);
				logMessage("Picked fault: " + pickedLine.getLabel());
				if(!addVertex(lineSurfaceVertex)) return true;
				model.win3d.win3dDisplay();
				return true;
			}
			else if(pickedEdge != null)
			{
				gridSectionPoint = new StsGridSectionPoint(pickedPoint, surface, true);
				StsSurfaceVertex newVertex = new StsSurfaceVertex(gridSectionPoint);
				newVertex.addEdgeAssociation(pickedEdge);
				if(!addVertex(newVertex)) return true;
                logMessage("Edge " + pickedEdge.getLabel() + " has been picked.");
				model.win3d.win3dDisplay();
				return true;
			}
		    else // nothing picked: make a point on the surface
			{
				if(gridPoint == null) return true;
				gridSectionPoint = new StsGridSectionPoint(gridPoint, true);
				if(nVerticesPicked == 0)
				{
					boolean pickOK = StsYesNoDialog.questionValue(model.win3d, "This is your first pick and it's not on an edge or fault-line.\n" +
                        "This means this will be a dying fault point.  Is this OK?");
                    if(!pickOK) return true;

					StsSurfaceVertex newVertex = new StsSurfaceVertex(gridSectionPoint);
					if(!addVertex(newVertex)) return true;
				}
				else
					addEdgePoint(gridSectionPoint);

				logMessage("Point added at: " + gridSectionPoint.getRowColLabel(surface));
				model.win3d.win3dDisplay();
				return true;
			}
		}
		catch(Exception e)
		{
			StsException.outputException("StsEdgeOnSurface.performMouseAction() failed.",
				e, StsException.WARNING);
			return false;
		}
	}

/*
   	public boolean performMouseAction(StsMouse mouse)
	{
        StsPickItem pickItem;
		StsSectionEdge pickedEdge;
		StsGridSectionPoint gridSectionPoint;
		StsMousePoint mousePoint;

		try
		{
			int leftButtonState = mouse.getButtonState(StsMouse.LEFT);

			if(leftButtonState == StsMouse.PRESSED || leftButtonState == StsMouse.DRAGGED)
			{
				mousePoint = mouse.getMousePoint();
				surface.getSurfacePosition(mousePoint, true);
				return true;
			}

			if(leftButtonState != StsMouse.RELEASED)  return true;

			// mouseButton RELEASED: process selection



			// haven't start edge yet: set action from Fault/Aux selection to "End Edge"
			if(edge == null) setActionButtonCommand("End Edge");

			mousePoint = mouse.getMousePoint();
			StsGridPoint gridPoint = surface.getSurfacePosition(mousePoint, true);
			if(gridPoint == null) return false;

			// see if a well is picked
			StsWell well = (StsWell)StsPick.pickClass3d(win3d, getConnectableWells(),
					StsPick.PICKSIZE_MEDIUM, StsPick.PICK_CLOSEST);

			if(well != null)
			{
				if(well.isFullyConnected())
				{
					StsMessage.printMessage("Sorry this fault is fully connected.\n" +
						"Select a point on an adjoining section or another fault.");
					return true;
				}

				StsSurfaceVertex wellVertex = well.getWellLine().getSectionEdgeVertex(surface);
				logMessage("Picked fault: " + well.getLabel());
				if(!addVertex(wellVertex)) return true;
				model.win3d.win3dDisplay();
				return true;
			}

			// no well picked: see if an edge is picked
			pickedEdge = (StsSectionEdge)StsPick.pickClass3d(win3d, getPickableEdges(),
						 StsPick.PICKSIZE_MEDIUM, StsPick.PICK_CLOSEST);

			if(pickedEdge != null)
			{
				logMessage("Intersected edge: " + pickedEdge.getLabel());
				pickItem = StsPick.pickItems[0];
				int nSegment = pickItem.names[1];
				StsPoint pickedPoint = pickedEdge.getPointOnLineNearestMouse(nSegment, mousePoint, win3d.glPanel3d);
				gridSectionPoint = new StsGridSectionPoint(pickedPoint, surface, true);
				StsSurfaceVertex newVertex = new StsSurfaceVertex(gridSectionPoint);
				newVertex.addEdgeAssociation(pickedEdge);
				if(!addVertex(newVertex)) return true;
				model.win3d.win3dDisplay();
				return true;
			}

			// nothing picked: make a point on the surface
			StsPoint point = new StsPoint(gridPoint.pxyz, 0.0f);
			gridSectionPoint = new StsGridSectionPoint(point, surface, true);
			if(nVerticesPicked == 0)
			{
				StsSurfaceVertex newVertex = new StsSurfaceVertex(gridSectionPoint);
				if(!addVertex(newVertex)) return true;
			}
			else
			    addEdgePoint(gridSectionPoint);

			logMessage("Point added at: " + gridSectionPoint.getRowColLabel(surface));
			model.win3d.win3dDisplay();
		    return true;
		}
		catch(Exception e)
		{
			StsException.outputException("StsEdgeOnSurface.performMouseAction() failed.",
				e, StsException.WARNING);
			return false;
		}
	}
*/
	private boolean addVertex(StsSurfaceVertex pickedVertex)
	{
	    if(!addVertexToEdge(pickedVertex)) return false;
		if(isEdgeComplete(pickedVertex)) return completeEdge();
		return true;
	}

	private boolean isEdgeComplete(StsSurfaceVertex pickedVertex)
	{
		// if we have picked the max number of vertices or we have picked more than 1
		// vertex and have now terminated on an edge, we are done with this edge
        return nVerticesPicked == nVerticesMax || nVerticesPicked > 1  && pickedVertex.isOnEdge();
//		return nVerticesPicked == nVerticesMax || nVerticesPicked > 1 && ( pickedVertex.isOnEdge() || pickedVertex.isOnLine() );
	}

	private boolean addVertexToEdge(StsSurfaceVertex pickedVertex)
	{
		if(edge == null)
		{
            edge = new StsSectionEdge(type, null, surface, nPointsPerVertex);
            edge.setDrawPoints(true);
            model.addDisplayableInstance(edge);
            sectionEdges.add(edge);
            edgeVertices = (StsSurfaceVertex[][])StsMath.arrayAddElement(edgeVertices, new StsSurfaceVertex[0]);
		}
		edgeVertices[nEdges] = (StsSurfaceVertex[])StsMath.arrayAddElement(edgeVertices[nEdges], pickedVertex);
		if(!checkEdgeDirection())
		{
		    new StsMessage(model.win3d, StsMessage.WARNING, "Sorry the vertex picked doesn't correspond to edges already picked.\n" +
						   "Deleting this edge. Try again or undo operation.");
			deleteEdge();
		    return false;
		}
		pickedVertex.setSurface(surface);
        pickedVertex.checkAddTimeOrDepth();
		nVerticesPicked++;
		edge.addEdgePoint(pickedVertex.getSurfacePoint());
		edge.setPointsFromEdgePoints();
		return true;
	}

	// pick is on the surface (not edge or fault);  add to new edge or existing edge
	private void addEdgePoint(StsGridSectionPoint edgePoint)
	{
		edge.addEdgePoint(edgePoint);
		edge.setPointsFromEdgePoints();
	}

    private boolean completeEdge()
    {
		if(edge == null) return false;

		StsGridSectionPoint edgePoint = (StsGridSectionPoint)edge.getEdgePointsList().getLast();
		StsSurfaceVertex lastVertex = edgePoint.getVertex();
        if(lastVertex == null)
        {
            boolean pickOK = StsYesNoDialog.questionValue(model.win3d, "Your last pick didn't terminate on an edge or fault-line.\n" +
                        "This means this will be a dying fault point.  Is this OK?");
            if(!pickOK) return false;
        }
		edge.setDrawPoints(false);
		if(!buildEdge()) return false;

		logMessage("  fault edge on " + surfaceName + " digitized successfully.");

		if(section == null) section = new StsSection(type, false);

		if(!addSectionEdge(edge))
		{
			new StsMessage(model.win3d, StsMessage.WARNING, "Failed to add this edge to section.\n" +
				"Deleting edge. Please try again.");
			deleteEdge();
			return false;
		}

		if(isSectionComplete()) completeSection();
		else setActionButtonCommand("End Section");

		clearEdge();
		return true;
    }

    private boolean buildEdge()
    {
		StsPoint firstEdgeVector, edgeVector;
		StsSurfaceVertex[] vertices;
		boolean ok;

		// if last pick is on surface, it has no vertex: add one

		StsGridSectionPoint edgePoint = (StsGridSectionPoint)edge.getEdgePointsList().getLast();
		StsSurfaceVertex lastVertex = edgePoint.getVertex();
		if(lastVertex == null)
		{
			lastVertex = new StsSurfaceVertex(edgePoint);
			if(!addVertexToEdge(lastVertex)) return false;
		}

		if(edgeDirection == StsParameters.NONE)
		{
			firstEdgeVector = getEdgeVector(0);
			edgeVector = getEdgeVector(nEdges);

			float dot = StsMath.dot(firstEdgeVector.v, edgeVector.v);
			if(dot < 0.0f) edgeDirection = StsParameters.MINUS;
			else edgeDirection = StsParameters.PLUS;
		}

		vertices = edgeVertices[nEdges];
		if(edgeDirection == StsParameters.MINUS)
		{
			StsMath.reverseOrder(vertices);
			edge.reverseEdgePoints();
		}

		edge.setPrevVertex(vertices[0]);
		edge.setNextVertex(vertices[vertices.length-1]);

		nEdges++;
		if(nEdges == 1) nVerticesMax = nVerticesPicked;

		surfacePicked[surfaceIndex] = true;
		return true;
    }

	private StsPoint getEdgeVector(int nEdge)
	{
		StsSurfaceVertex[] vertices = edgeVertices[nEdge];
		StsPoint point0 = vertices[0].getPoint();
		StsPoint point1 = vertices[vertices.length-1].getPoint();
		StsPoint point = StsPoint.subPointsStatic(point1, point0);
//		point.subPoints(point1, point0);
		return point;
	}

	private boolean checkEdgeDirection()
	{
		int nEdges = edgeVertices.length;
		if(nEdges < 2) return true;
		StsObject[] firstEdgeAssocs = getVertexAssociations(edgeVertices[0]);
		if(firstEdgeAssocs == null) return true;
		StsObject[] currentEdgeAssocs = getVertexAssociations(edgeVertices[nEdges-1]);
		if(currentEdgeAssocs == null) return true;

		int nVerticesPicked = currentEdgeAssocs.length;
		int nVerticesMax = firstEdgeAssocs.length;

		// check forward direction
		boolean forwardMatches = true;
		for(int n = 0; n < nVerticesPicked; n++)
		{
			if(firstEdgeAssocs[n] != currentEdgeAssocs[n])
			{
				forwardMatches = false;
				break;
			}
		}
		boolean reverseMatches = true;
		for(int n = 0; n < nVerticesPicked; n++)
		{
			if(firstEdgeAssocs[nVerticesMax-1-n] != currentEdgeAssocs[n])
			{
				reverseMatches = false;
				break;
			}
		}

		if(!forwardMatches && !reverseMatches) return false;

		if(forwardMatches && !reverseMatches) edgeDirection = StsParameters.PLUS;
		else if(reverseMatches && !forwardMatches) edgeDirection = StsParameters.MINUS;
		return true;
	}

	private StsObject[] getVertexAssociations(StsSurfaceVertex[] vertices)
	{
		if(vertices == null) return null;
		int nVertices = vertices.length;
		StsObject[] associations = new StsObject[nVertices];
		for(int n = 0; n < nVertices; n++)
			associations[n] = vertices[n].getAssociatedSectionOrLine();
		return associations;
	}
/*
	private void checkReverseInitialEdgeOrder()
	{
		if(edge == null) return;
		if(edge.checkReverseEdgeOrder()) StsMath.reverseOrder(edgeVertices[0]);
	}
*/
/*
    public boolean checkReverseEdgeOrder()
    {
		StsSectionEdge otherEdge;
		StsSection otherSection;

		boolean reversedEdge = false;
        // connection order is OK, NOT_OK (reversed), or NONE (no other connection)
        int minusEndOrder = edge.getConnectedEdgeOrder(MINUS);
        int plusEndOrder = edge.getConnectedEdgeOrder(PLUS);

        if(minusEndOrder == NOT_OK)
        {
            if(plusEndOrder == OK)
			{
			    otherEdge = edge.getConnectedEdge(MINUS);
				otherSection = otherEdge.getSection();
				if(otherSection != null) otherSection.reverse();
			}
            else
			{
                reverseEdge();
				reversedEdge = true;
			}
        }
        else if(plusEndOrder == NOT_OK)
        {
            if(minusEndOrder == OK)
			{
			    otherEdge = edge.getConnectedEdge(PLUS);
				otherSection = otherEdge.getSection();
				if(otherSection != null) otherSection.reverse();
			}
            else
			{
                reverseEdge();
				reversedEdge = true;
			}
        }
		return reversedEdge;
    }
*/
	// mainDebug purposes only: move intermediate edge points to nearest grid row-col
	// with only one edge, this means that section col lines will also fall on grid row & col lines
	// which results in grid row & col edges emanating from same edge  point
	// which is what we want to test
	private void edgeAdjustToGridPoints()
	{
		StsGridSectionPoint edgePoint;
		float rowF, colF;

		StsList edgePoints = edge.getEdgePointsList();
		int nEdgePoints = edgePoints.getSize();
		for(int n = 1; n < nEdgePoints-1; n++)
		{
			edgePoint = (StsGridSectionPoint)edgePoints.getElement(n);
			int row = StsMath.roundOffInteger(edgePoint.getGridRowF());
			int col = StsMath.roundOffInteger(edgePoint.getGridColF());
			edgePoint.setRowOrColIndex(null, StsParameters.ROW, row);
			edgePoint.setRowOrColIndex(null, StsParameters.COL, col);
		}
		edge.setPointsFromEdgePoints();
	}
/*
  	private boolean completeEdgeAndSection()
    {
		try
		{
			boolean ok = true;

//			status.textOnly();
			surface.setIsVisible(visible);
			model.viewProperties.restore();

			if(edge != null)
			{
				edge.setDrawPoints(false);
				if(!completeEdge())
				{
					edge.delete();
					logMessage(" fault edge on " + surfaceName + " not digitized successfully.", 1000, true);
					ok = false;
				}
				else
				{
					model.setModelState(StsModelState.SECTIONS);
					logMessage("  fault edge on " + surfaceName + " digitized successfully.", 1000, true);
					ok = true;
				}

				if(section == null) section = new StsSection(model, type);

				if(!addSectionEdge(edge))
					logMessage("Failed to add this edge to section.");

				incrementSurface();
			}

			if(section != null)
			{
				if(completed)
					completeSection();
				else
					setActionButtonCommand("End Section");
			}
			return ok;
		}
		catch(Exception e)
		{
			StsException.outputException("StsEdgeOnSurface.completeEdgeAndSection() failed.",
				e, StsException.WARNING);
			return false;
		}
		finally
		{
			clearEdge();
		}
    }
*/
	private boolean addSectionEdge(StsSectionEdge edge)
	{
		StsObjectRefList sectionEdges = section.getSectionEdges();
		int nSectionEdges = sectionEdges.getSize();
		if(nSectionEdges == 0)
		    return section.addSectionEdge(edge);
		else
		{
			int surfaceIndex1 = StsParameters.NO_MATCH; // -1
			for(int n = 0; n < nSectionEdges; n++)
			{
				StsSectionEdge sectionEdge = (StsSectionEdge)sectionEdges.getElement(n);
				int surfaceIndex0 = surfaceIndex1;
				surfaceIndex1 = getSurfaceIndex(sectionEdge.getSurface());
				if(surfaceIndex > surfaceIndex0 && surfaceIndex < surfaceIndex1)
				{
				    sectionEdges.insertBefore(n, edge);
					return true;
				}
			}
			return section.addSectionEdge(edge);
		}
	}

	private int getSurfaceIndex(StsXYSurfaceGridable surface)
	{
		if(surfaces == null) return StsParameters.NO_MATCH;
        for(int n = 0; n < surfaces.length; n++)
            if(surface == surfaces[n]) return n;
        return StsParameters.NO_MATCH;
	}

	private void clearEdge()
	{
		edge = null;
		nVerticesPicked = 0;
		edgeDirection = StsParameters.NONE;
	}

	private boolean completeSection()
	{
		try
		{
			// would prefer delaying this addToModel so that we don't have 2
			// db updates for section in this method;
			// however, we need the section.index set so that this section
			// can be added to fault.connectedSections objRefList which
			// require that obj being added has a persistent index (>= 0).
//			section.addToModel();

			StsFaultLine[] faults = addSectionFaults(section);
			if(faults == null)
            {
                deleteSection();
                return false;
            }
			if(!addSectionToIntermediateVertices(section)) deleteSection();
			adjustConnectedSections(section);
			section.addToModel();
			addFaultsToSection(section, faults);
			section.completeSection();
			// section.constructSection();
			model.refreshObjectPanel(); // rebuilds object tree
			initializeSectionPick();
			repeatAction();
			this.nSectionsBuilt++;
            model.getProject().setModelZDomainToCurrent();

            return true;
		}
		catch(Exception e)
		{
			StsException.outputException("StsEdgeOnSurface.completeSection() failed.",
				e, StsException.WARNING);
			deleteSection();
			return false;
		}
		finally
		{
			setButtonEdgeTypeCommand();
		}
	}

	private void addFaultsToSection(StsSection section, StsFaultLine[] faults)
	{
		section.addLines(faults);
//		section.dbFieldChanged("lines", section.getLines());
	}

    private void deleteSection()
	{
		if(section != null)
		{
			StsObjectRefList.deleteAll(section.getLines());
		    section.delete();
            removeDisplayableEdges();
		}
		initializeSectionPick();
		model.win3d.win3dDisplay();

	}

	private StsFaultLine[] addSectionFaults(StsSection section)
	{
		StsFaultLine[] faults = new StsFaultLine[2];
		StsSurfaceVertex[] vertices;
		StsSection onSection;

		try
	    {
			StsSeismicVelocityModel velocityModel = model.getProject().getSeismicVelocityModel();
			vertices = getVerticesAtIndex(0);
            if(vertices == null || vertices[0] == null) return null;
			onSection = vertices[0].onSection(); // we are assuming all vertices on same section
			faults[0] = getCommonVertexFault(vertices);
			if(faults[0] == null) faults[0] = constructFaultFromVertices(vertices);
			if(faults[0] == null) return null;
//			if(velocityModel != null) faults[0].adjustTimeOrDepth(velocityModel);
//			section.addLine(fault);

			vertices = getVerticesAtIndex(nVerticesMax-1);
			onSection = vertices[0].onSection();
			faults[1] = getCommonVertexFault(vertices);
			if(faults[1] == null) faults[1] = constructFaultFromVertices(vertices);
			if(faults[1] == null) return null;
//			if(velocityModel != null) faults[1].adjustTimeOrDepth(velocityModel);
//			section.addLine(fault);
			return faults;
	    }
	    catch(Exception e)
		{
			StsException.outputException("StsEdgeOnSurface.getConstructSectionFaults() failed.",
				e, StsException.WARNING);
			return null;
		}
	}

	private StsFaultLine constructFaultFromVertices(StsSurfaceVertex[] vertices)
	{
		try
		{
			if(vertices == null) return null;
			int nVertices = vertices.length;
		    if(nVertices <= 0) return null;

//            StsGridSectionPoint[] gridPoints = new StsGridSectionPoint[nVertices];
//            for(int n = 0; n < nVertices; n++)
//                gridPoints[n] = vertices[n].getSurfacePoint();
			StsFaultLine line = StsFaultLine.buildFault(vertices);
            line.extendEnds(); // this constructs the line including the extended ends
			StsSection onSection = vertices[0].getAssociatedSection();
			if(onSection != null)
            {
                line.setOnSection(onSection);
                line.projectToSection();
            }
            line.initialized = true;
//            line.checkAdjustDepthOrTimePoints();
			return line;
		}
		catch(Exception e)
		{
			StsException.outputException("StsEdgeOnSurface.getConstructEndFault() failed.",
				e, StsException.WARNING);
			return null;
		}
	}

	private StsSurfaceVertex[] getVerticesAtIndex(int vertexIndex)
	{
		try
		{
            if(edgeVertices == null || edgeVertices.length < nSurfaces) return null;

			StsSurfaceVertex[] vertices = new StsSurfaceVertex[0];
			for(int s = 0; s < nSurfaces; s++)
			{
				if(edgeVertices[s] != null && edgeVertices[s].length > vertexIndex && edgeVertices[s][vertexIndex] != null)
					vertices = (StsSurfaceVertex[])StsMath.arrayAddElement(vertices, edgeVertices[s][vertexIndex]);
			}
			return vertices;
		}
		catch(Exception e)
		{
			StsException.outputException("StsEdgeOnSurface.getVerticesAtIndex() failed.",
				e, StsException.WARNING);
			return new StsSurfaceVertex[0];
		}
	}

	private StsFaultLine getCommonVertexFault(StsSurfaceVertex[] vertices)
	{
		StsFaultLine fault = null;
		int nVertices = vertices.length;
		for(int n = 0; n < nVertices; n++)
		{
			StsLine vertexLine = vertices[n].getSectionLine();
            if(vertexLine != null && vertexLine instanceof StsFaultLine)
            {
                StsFaultLine vertexFault = (StsFaultLine)vertexLine;
				if(fault != null && fault != vertexFault)
				{
					StsException.systemError("StsEdgeOnSurface.getCommonVertexFault() failed." +
						" Found different faults at this vertex location: " +
						fault.getLabel() + " and " + vertexFault.getLabel());
					return null;
				}
				fault = vertexFault;
			}
		}
		return fault;
	}

	private boolean addSectionToIntermediateVertices(StsSection section)
	{
		try
		{
			if(section == null) return false;

			for(int n = 1; n < nVerticesPicked-1; n++)
			{
				StsSurfaceVertex[] vertices = getVerticesAtIndex(n);
				if(vertices == null)
				{
					StsException.systemError("StsEdgeOnSurface.addSectionToIntermediateVertices() failed." +
						" Couldn't find intermediate vertices at index: " + n);
					return false;
				}
				StsLine fault = getCommonVertexFault(vertices);
				if(fault == null)
				{
					StsException.systemError("StsEdgeOnSurface.addSectionToIntermediateVertices() failed." +
						" Couldn't find common fault for vertices at index: " + n);
					return false;
				}
				for(int v = 0; v < vertices.length; v++)
					vertices[v].setSectionLine(fault);

				if(!fault.addOnSection(section)) return false; // this existing fault is attached to new section
			}
			return true;
		}
		catch(Exception e)
		{
			StsException.outputException("StsEdgeOnSurface.getConstructEndFault() failed.",
				e, StsException.WARNING);
			return false;
		}
	}

	private void adjustConnectedSections(StsSection section)
	{
		if(!section.hasLines()) return;
		StsClass sections = model.getCreateStsClass(StsSection.class);
		int nSections = sections.getSize();
		for(int n = 0; n < nSections; n++)
			((StsSection)sections.getElement(n)).initialized = true;

		ConnectedSection rootSection = new ConnectedSection(section);
		section.initialized = section.isEndSectionsInitialized();

		boolean allInitialized = false;
		int maxIters = 3;
		int iter = 0;
		while(!allInitialized && iter++ < maxIters)
		    allInitialized = rootSection.reinitialize();
	}

	class ConnectedSection
	{
		StsSection section;
		StsLine[] lines = new StsLine[0];
		ConnectedSection[] connectedSections = new ConnectedSection[0];

		ConnectedSection(StsSection section)
		{
			this.section = section;
			section.initialized = false;
			addConnectedSideSections(section.getLeftLines());
			addConnectedSideSections(section.getRightLines());
		}

		void addConnectedSideSections(StsObjectRefList sectionLines)
		{
			if(sectionLines == null) return;
			int nLines = sectionLines.getSize();
			for(int n = 0; n < nLines; n++)
			{
				StsLine line = (StsLine)sectionLines.getElement(n);
				lines = (StsLine[])StsMath.arrayAddElement(lines, line);

				StsSection sideSection = line.getOnlyConnectedSection();
				if(sideSection != null && sideSection.initialized)
					connectedSections = (ConnectedSection[])StsMath.arrayAddElement(connectedSections,
										    new ConnectedSection(sideSection));
			}
		}

		boolean reinitialize()
		{
			boolean initialized = true;

			if(!section.initialized)
			{
				initialized = initialized && section.reinitialize();
				reinitializeLines();
			}

			int nConnectedSections = connectedSections.length;
			for(int n = 0; n < nConnectedSections; n++)
				initialized = initialized && connectedSections[n].reinitialize();

			return initialized;
		}

		void reinitializeLines()
		{
			int nLines = lines.length;
			for(int n = 0; n < nLines; n++)
				lines[n].reinitialize();
		}
	}



    private StsObject[] getPickableEdges()
    {
        return getPickableEdgesExcept(edge);
    }

    private StsObject[] getPickableEdgesExcept(StsSectionEdge exceptEdge)
    {
        StsClass sectionEdges = model.getCreateStsClass(StsSectionEdge.class);
        int nSectionEdges = sectionEdges.getSize();
        StsObject[] selectedEdges = new StsObject[nSectionEdges];
        int nSelectedEdges = 0;
        for(int n = 0; n < nSectionEdges; n++)
        {
            StsSectionEdge edge = (StsSectionEdge)sectionEdges.getElement(n);

//            if(edge == exceptEdge) //mainDebug
//                StsException.systemDebug("Found edge matching vertexEdge: " + edge.index());

            if(edge.getType() != StsParameters.REFERENCE && edge != exceptEdge && edge.getSurface() == surface)
                selectedEdges[nSelectedEdges++] = edge;
        }
        return selectedEdges;
    }

	// connectable lines here must be: displayed, faults/pseudos only, not fully-connected
    private StsFaultLine[] getConnectableLines()
    {
        StsMainObject[] faultLines = model.getVisibleObjectList(StsFaultLine.class);
		int nFaultLines = faultLines.length;
		StsFaultLine[] connectableLines = new StsFaultLine[nFaultLines];
		int nConnectableLines = 0;
		for(int n = nFaultLines-1; n >= 0; n--)
		{
            StsFaultLine faultLine = (StsFaultLine)faultLines[n];
			if(!faultLine.isFullyConnected()) connectableLines[nConnectableLines++] = faultLine;
		}
		return (StsFaultLine[])StsMath.trimArray(connectableLines, nConnectableLines);
    }

	private void setButtonEdgeTypeCommand()
	{
		if(type == StsParameters.FAULT)
			setActionButtonCommand("Fault");
		else
			setActionButtonCommand("Auxiliary");
	}

	public boolean end()
    {
        boolean ok = true;
		setIsRepeatable(false);
        statusArea.textOnly();
		surface.setIsVisible(visible);
//        model.viewProperties.restore();
		if(section != null)
		{
			if (edge != null) if(!completeEdge()) ok = false;
			if(!completeSection()) ok = false;
		}
		else
			ok = true;
//			ok = completeSection();

        model.setActionStatus(StsBuildFrame.class.getName(), StsModel.STARTED);
        clearToolbar();
        removeDisplayableEdges();
        model.win3d.win3dDisplay();
        return ok;
    }

    private void removeDisplayableEdges()
    {
        int nEdges = sectionEdges.getSize();
        for(int n = 0; n < nEdges; n++)
        {
            StsSectionEdge edge = (StsSectionEdge)sectionEdges.getElement(n);
            model.removeDisplayableInstance(edge);
        }
    }

	// increment surface; if all surfaces picked, return true (completed), otherwise false
	private boolean isSectionComplete()
	{
		for(int n = 0; n < nSurfaces; n++)
		{
			surfaceIndex = (surfaceIndex+1)%nSurfaces;
			if(!surfacePicked[surfaceIndex])
			{
				setSurface();
				return false;
			}
		}
		return true;
	}
}


