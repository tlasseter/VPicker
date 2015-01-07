package com.Sts.Actions.Wizards.Fracture;

import com.Sts.Actions.Build.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.awt.event.*;

public class StsBuildFractureSet extends StsWizardStep
{
    StsBuildFractureSetPanel panel;
    StsHeaderPanel header;
    
    static public final int NONE = 0;    
    static public final int TYPE_EDGE = 1;
    static public final int TYPE_VERTEX = 2;
    static public int keyType = NONE;
       
    static public final int PICKING = 1;
    static public final int DELETING = 2;    
    static public int mode = PICKING;
    
    private StsFractureSet currentFractureSet;
    private StsCursorZSurface cursorZSurface; // current surface
	private StsCursorZSurface[] cursorZSurfaces = new StsCursorZSurface[0]; // array of top-down surfaces
	private int nCursorZs; // number of surfaces
//	private boolean[] zCursorPicked; // indicates surface already picked
	private int nEdges = 0;  // number of edges constructed
	private StsSurfaceVertex[][] edgeVertices; // [nEdges][nVertices] array
	private int cursorZIndex = 0; // curent index of edge
	private int nVerticesPicked = 0; // number of vertices picked for this edge
	private int firstVertexIndex = -1;  // vertex index in edgeVertices of first pick
	private int nVerticesMax = 0; // max nVertices determined from first edge picked
    private boolean visible = true;
    private StsSectionEdge edge;
    private StsObjectList sectionEdges = new StsObjectList(2, 2);
    private StsFracture section;
	private int nSectionsBuilt = 0;
    private boolean ok = true;
	private int edgeDirection; // NONE, PLUS, or MINUS
    private boolean addDepth = false;

    private byte type = StsParameters.FAULT;  // FAULT, or AUXILIARY
	static private int nPointsPerVertex = 1;

   	/** create a new edge on surface and attach a mouse listener to it */
 	public StsBuildFractureSet(StsWizard wizard)
    {
    	super(wizard);
        panel = new StsBuildFractureSetPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Build Fracture Set");
        header.setSubtitle("Defining Fracture Set in 3D");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Fracture");        
        header.setInfoText(wizardDialog,"(1) Set the Z Slice to the desired location.\n" +
        		"(2) Digitize the edge on the slice.\n" +
        		"(3) Press the End Edge Button when the edge is complete.\n" +
        		" ***** Repeat steps 1-3 until all edges are defined\n" +
        		"(4) Press the End Section Button to complete the current fracture.\n" +
        		" ***** Repeat steps 1-4 until all sections are complete\n" +
        		" ***** The Delete Edge Button can be used to delete the current edge only\n"+
        		"(5) Press the Finish>> Button when the fracture set is complete.\n");
        header.setLink("http://www.s2ssystems.com/marketing/s2ssystems/BuildFracture.html");
    }

	public boolean start()
    {
		if( !actionOK() )  return false;

        try
        {
            statusArea.setTitle("Pick fault edges on Z cursor: ");
            currentFractureSet = (StsFractureSet)((StsFractureWizard)wizard).getSelectedFractureSet();
			cursorZSurfaces = currentFractureSet.getCursorZSurfaces();
            initializeSectionPick();
            turnOffFaultLines();
			logMessage("Select a fault edge on surface or 'End Edge' button to skip.");

            StsProject project = model.getProject();
            byte zDomainSupported = project.getZDomainSupported();
            addDepth = zDomainSupported == StsProject.TD_TIME_DEPTH || zDomainSupported == StsProject.TD_APPROX_DEPTH_AND_DEPTH;
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

    private void turnOffFaultLines()
    {
        StsFaultLineClass faultLineClass = (StsFaultLineClass)model.getStsClass(StsFaultLine.class);
        faultLineClass.setDisplayLines(false);
    }

    private void initializeSectionPick()
	{
		nCursorZs = cursorZSurfaces.length;
//		zCursorPicked = new boolean[nCursorZs]; // assume initialized to false
		edgeVertices = new StsSurfaceVertex[0][];
		section = null;
	    nEdges = 0;
	    firstVertexIndex = -1;
	    nVerticesMax = 0;
		cursorZIndex = 0;
        setCursorZ();
        panel.endEdgeButton.setEnabled(false);
        panel.deleteCurrentEdgeButton.setEnabled(false);
        panel.endSectionButton.setEnabled(false);
    }

	private void setCursorZ()
	{
        if(cursorZSurface != null)
            cursorZSurface.toggleSurfacePickingOff();
        if(nCursorZs == 0) return;
        setCursorZSurface();
        model.win3d.adjustCursorAndSlider(StsCursor3d.ZDIR, cursorZSurface.getZ());
    }

    private void setCursorZSurface()
    {
        cursorZSurface = cursorZSurfaces[cursorZIndex];
        cursorZSurface.toggleSurfacePickingOn();
        cursorZSurface.setIsVisible(true);
    }

	private boolean actionOK()
    {
	    return true;
    }

    public void endEdge()
    {
        if(!completeEdge()) 
        	return;
        panel.endEdgeButton.setEnabled(false);
        panel.deleteCurrentEdgeButton.setEnabled(false);
        unlockZCursor();
        if(nEdges > 1) panel.endSectionButton.setEnabled(true);
        model.win3d.win3dDisplay();
    }

    private void unlockZCursor()
    {
        model.win3d.cursor3dPanel.unlockSlider(StsCursor3d.ZDIR);
    }

    private void lockZCursor()
    {
        model.win3d.cursor3dPanel.lockSlider(StsCursor3d.ZDIR);
    }

    public boolean endSection()
    {
        if(!completeSection()) return false;
        panel.endSectionButton.setEnabled(false);
        model.win3d.win3dDisplay();
        return true;
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
		model.win3d.win3dDisplay();
	}

    // mouse action for 3d window
   	public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel)
	{
		StsFaultLine pickedLine;
		StsSectionEdge pickedEdge;
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
                StsPoint cursorZPoint = model.win3d.getCursor3d().getPointInCursorPlane(glPanel3d, StsCursor3d.ZDIR, mouse);
				return true;
			}

			if(leftButtonState != StsMouse.RELEASED)  return true;

			// mouseButton RELEASED: process selection

			// haven't start edge yet: set action from Fault/Aux selection to "End Edge"
//			panel.endEdgeButton.setEnabled(false);
            // get intersection point on surface (if it exists)
            StsPoint cursorZPoint = model.win3d.getCursor3d().getPointInCursorPlane(glPanel3d, StsCursor3d.ZDIR, mouse);
            checkCursorZSurface();
            cursorZPoint = cursorZSurface.getXYZTPoint(cursorZPoint);
            StsGridPoint gridPoint = new StsGridPoint(cursorZPoint, model.getProject().getGridDefinition());
//            if(gridPoint == null) return true;

			// see if a line is picked
			pickedLine = (StsFaultLine)StsJOGLPick.pickClass3d(glPanel, (StsObject[])getConnectableLines(), StsJOGLPick.PICKSIZE_MEDIUM, StsJOGLPick.PICK_CLOSEST);
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
					pickedPoint = pickedLine.getPointOnLineNearestMouse(nSegment, mouse, (StsGLPanel3d)glPanel);

                    if(gridPoint != null)
                    {
                        float[] xyz = pickedPoint.getXYZorT();
					    screenZLine = glPanel.getScreenZ(xyz);
//			            float screenZSurface = glPanel.getScreenZ(gridPoint.getXYZorT());
//				        if(screenZLine >= screenZSurface) pickedLine = null;
                    }
				}
			}

			pickedEdge = (StsSectionEdge)StsJOGLPick.pickClass3d( glPanel, getPickableEdges(), StsJOGLPick.PICKSIZE_MEDIUM, StsJOGLPick.PICK_CLOSEST);
			if(pickedEdge != null)
			{
				pickItem = StsJOGLPick.pickItems[0];
				int nSegment = pickItem.names[1];
				pickedPoint = pickedEdge.getPointOnLineNearestMouse(nSegment, mouse, (StsGLPanel3d)glPanel);
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
				StsSurfaceVertex lineSurfaceVertex = pickedLine.getSectionEdgeVertex(cursorZSurface);
				logMessage("Picked fault: " + pickedLine.getLabel());
				if(!addVertex(lineSurfaceVertex)) return true;
				model.win3d.win3dDisplay();
				return true;
			}
			else if(pickedEdge != null)
			{
				gridSectionPoint = new StsGridSectionPoint(pickedPoint, cursorZSurface, true);
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
            /*
                    boolean pickOK = StsMessage.questionValue(model.win3d, "This is your first pick and it's not on an edge or fault-line.\n" +
                        "This means this will be a dying fault point.  Is this OK?");
                    if(!pickOK) return true;
             */
					StsSurfaceVertex newVertex = new StsSurfaceVertex(gridSectionPoint);
					if(!addVertex(newVertex)) 
					{
						panel.deleteCurrentEdgeButton.setEnabled(true);						
						return true;
					}
				}
				else
				{
					addEdgePoint(gridSectionPoint);
					panel.deleteCurrentEdgeButton.setEnabled(true);
				}
            	
                int nEdgePoints = edge.getPoints().length;
                if(nEdgePoints == 1)
                {
                	lockZCursor();
                }
                else 
                {
                	if(nEdgePoints == 2) 
                		panel.endEdgeButton.setEnabled(true);
                }
                logMessage("Point added at: " + gridSectionPoint.getRowColLabel(cursorZSurface));
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

    private void checkCursorZSurface()
    {
        float currentZ = model.getCursor3d().getCurrentDirCoordinate(StsCursor3d.ZDIR);
        for(int n = 0; n < nCursorZs; n++)
        {
            float cursorZ = cursorZSurfaces[n].getZ();
            if(StsMath.sameAs(currentZ, cursorZ))
            {
                cursorZSurface = cursorZSurfaces[n];
                cursorZIndex = n;
                return;
            }
            if(currentZ < cursorZ)
            {
                cursorZSurface = new StsCursorZSurface(currentZ);
                cursorZSurfaces = (StsCursorZSurface[])StsMath.arrayInsertElementBefore(cursorZSurfaces, cursorZSurface, n);
                cursorZIndex = n;
                nCursorZs++;
                return;
            }
        }
        cursorZSurface = new StsCursorZSurface(currentZ);
        cursorZIndex = cursorZSurfaces.length;
        cursorZSurfaces = (StsCursorZSurface[])StsMath.arrayInsertElementBefore(cursorZSurfaces, cursorZSurface, cursorZIndex);
        nCursorZs++;
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

			mouse.clearButtonState(StsMouse.LEFT);

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
            edge = new StsSectionEdge(type, null, cursorZSurface, nPointsPerVertex);
            edge.setDrawPoints(true);
            model.addDisplayableInstance(edge);
            sectionEdges.add(edge);
            edgeVertices = (StsSurfaceVertex[][])StsMath.arrayAddElement(edgeVertices, new StsSurfaceVertex[0]);
		}
		edgeVertices[nEdges] = (StsSurfaceVertex[])StsMath.arrayAddElement(edgeVertices[nEdges], pickedVertex);
		if(!checkEdgeDirection())
		{
		    new StsMessage(model.win3d, StsMessage.WARNING, "Sorry the vertex picked doesn't correspond to edges already picked.\n" +
						   "Must be picked in the same direction as other edges on section" +
						   "Deleting this edge. Try again or undo operation.");
			deleteEdge();
		    return false;
		}
		pickedVertex.setSurface(cursorZSurface);
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
    /*
        if(lastVertex == null)
        {
            boolean pickOK = StsMessage.questionValue(model.win3d, "Your last pick didn't terminate on an edge or fault-line.\n" +
                        "This means this will be a dying fault point.  Is this OK?");
            if(!pickOK) return false;
        }
    */
		edge.setDrawPoints(false);
		if(!buildEdge()) return false;

		logMessage("  fault edge on cursor Z surface at z" + cursorZSurface.getZ() + " digitized successfully.");

		if(section == null) 
		{
			section = new StsFracture(currentFractureSet);
			section.setStsColor((StsColor)((StsFractureWizard)wizard).getFractureSetColor());
		}

		if(!addSectionEdge(edge))
		{
			new StsMessage(model.win3d, StsMessage.WARNING, "Failed to add this edge to section.\n" +
				"Deleting edge. Please try again.");
			deleteEdge();
			return false;
		}
		// Add cursorZSurface reference to fractureSet
		currentFractureSet.addCursorZSurface(cursorZSurface);
		if(isSectionComplete()) endSection();
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

//		zCursorPicked[cursorZIndex] = true;
		return true;
    }

	private StsPoint getEdgeVector(int nEdge)
	{
		StsSurfaceVertex[] vertices = edgeVertices[nEdge];
		StsPoint point0 = vertices[0].getPoint();
		StsPoint point1 = vertices[vertices.length-1].getPoint();
		StsPoint point = StsPoint.subPointsStatic(point1, point0);
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
				if(cursorZIndex > surfaceIndex0 && cursorZIndex < surfaceIndex1)
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
		if(cursorZSurfaces == null) return StsParameters.NO_MATCH;
        for(int n = 0; n < cursorZSurfaces.length; n++)
            if(surface == cursorZSurfaces[n]) return n;
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
			section.addToSet();
			model.refreshObjectPanel(); // rebuilds object tree
			initializeSectionPick();
			repeatAction();
			this.nSectionsBuilt++;
            wizard.enableFinish();
			return true;
		}
		catch(Exception e)
		{
			StsException.outputException("StsEdgeOnSurface.completeSection() failed.",
				e, StsException.WARNING);
			deleteSection();
			return false;
		}
	}

	private void addFaultsToSection(StsSection section, StsFaultLine[] faults)
	{
		section.addLines(faults);
//		section.dbFieldChanged("lines", section.getLines());
	}

    public void deleteSection()
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
			if(velocityModel != null) faults[0].adjustTimeOrDepth(velocityModel);
//			section.addLine(fault);

			vertices = getVerticesAtIndex(nVerticesMax-1);
			onSection = vertices[0].onSection();
			faults[1] = getCommonVertexFault(vertices);
			if(faults[1] == null) faults[1] = constructFaultFromVertices(vertices);
			if(faults[1] == null) return null;
			if(velocityModel != null) faults[1].adjustTimeOrDepth(velocityModel);
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
            if(edgeVertices == null || edgeVertices.length < nCursorZs) return null;

			StsSurfaceVertex[] vertices = new StsSurfaceVertex[0];
			for(int s = 0; s < nCursorZs; s++)
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

		StsBuildFractureSet.ConnectedSection rootSection = new StsBuildFractureSet.ConnectedSection(section);
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
		StsBuildFractureSet.ConnectedSection[] connectedSections = new StsBuildFractureSet.ConnectedSection[0];

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
					connectedSections = (StsBuildFractureSet.ConnectedSection[])StsMath.arrayAddElement(connectedSections,
										    new StsBuildFractureSet.ConnectedSection(sideSection));
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

            if(edge.getType() != StsParameters.REFERENCE && edge != exceptEdge && edge.getSurface() == cursorZSurface)
                selectedEdges[nSelectedEdges++] = edge;
        }
        return selectedEdges;
    }

	// connectable lines here must be: isVisible, faults/pseudos only, not fully-connected
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

	public boolean end()
    {
        boolean ok = true;
		setIsRepeatable(false);
        statusArea.textOnly();
		cursorZSurface.setIsVisible(visible);
//        model.viewProperties.restore();
		if(section != null)
		{
			if (edge != null) if(!completeEdge()) ok = false;
			if(!endSection()) ok = false;
		}
		else
			ok = true;
//			ok = completeSection();

        model.setActionStatus(StsBuildFrame.class.getName(), StsModel.STARTED);
        clearToolbar();
        removeDisplayableEdges();
        model.viewObjectChangedAndRepaint(this, currentFractureSet);
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
        if(cursorZIndex == nCursorZs-1) return false;
        cursorZIndex++;
	    setCursorZ();
		return false;
	}
	
    public boolean keyReleased(KeyEvent e, StsMouse mouse, StsGLPanel glPanel)
    {
        if(mouse.isButtonDown()) return false;
        char keyReleased = e.getKeyChar();
        if(keyReleased == 'D' || keyReleased == 'd') // delete the current vertex or edge
        {
            if(StsYesNoDialog.questionValue(wizard.frame, "Are you sure you want to delete the current edge?"))
            {

            }
            return true;
        }
        return false;
    }
    
    public boolean deleteCurrentEdge()
    {
    	deleteEdge();
    	unlockZCursor();
    	StsMessageFiles.infoMessage("Edge successfully removed, reset to previous edge."); 
    	return true;
    }
    
    public boolean setCurrentEdgeToLast()
    {
    	edge = (StsSectionEdge)sectionEdges.getLast(); 
    	return true;
    }    
}
