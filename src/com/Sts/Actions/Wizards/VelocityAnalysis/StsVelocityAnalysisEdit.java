package com.Sts.Actions.Wizards.VelocityAnalysis;

import com.Sts.Actions.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Reflect.*;
import com.Sts.Types.*;
import com.Sts.Types.PreStack.*;
import com.Sts.UI.*;
import com.Sts.UI.Toolbars.*;
import com.Sts.Utilities.*;

import java.awt.event.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version 1.1
 */

public class StsVelocityAnalysisEdit extends StsAction
{
	StsPreStackLineSetClass lineSetClass;
	StsPreStackLineSet lineSet;
	StsPreStackVelocityModel velocityVolume;
	StsVelocityProfile velocityProfile;
    StsPreStackLine2d currentLine2d = null;
    StsPoint vertex;
	StsPoint pick;
    int typePicked = TYPE_NONE;
	int indexPicked = StsGather.NONE_INDEX;
	int nPickedCurve = -1;
//    boolean vertexChanged = false;
    static public final int TYPE_NONE = 0;
	static public final int TYPE_GATHER_EDGE = 1;
	static public final int TYPE_GATHER_VERTEX = 2;


    static public final byte CHANGE_NONE = StsVelocityProfile.CHANGE_NONE;
    static public final byte CHANGE_POINT = StsVelocityProfile.CHANGE_POINT;
    static public final byte CHANGE_MUTE = StsVelocityProfile.CHANGE_MUTE;

    static final boolean debug = false;

	public StsVelocityAnalysisEdit(StsActionManager actionManager)
	{
		super(actionManager);
	}

	public boolean start()
	{
		try
		{
			lineSet = StsPreStackLineSetClass.currentProjectPreStackLineSet;
            if(lineSet == null)return false;
 			velocityVolume = lineSet.getVelocityModel();
			if(velocityVolume == null)
			{
				new StsMessage(model.win3d, StsMessage.WARNING, "Failed to find or construct an initial velocity model.");
				return false;
			}
            lineSetClass = (StsPreStackLineSetClass)lineSet.getStsClass();
			lineSet.setMode(StsPreStackLineSet.EDIT_MODE);
            model.viewObjectRepaint(this, lineSet.getSemblanceDisplayProperties());
			logMessage("Pick velocity profile on semblanceBytes or gather displays.");
			return true;
		}
		catch(Exception e)
		{
			StsException.outputException("StsVelocityAnalysisEdit.start() failed.",
										 e, StsException.WARNING);
			return false;
		}
	}

    private void checkSetCurrentLine2d(StsVelocityProfile velocityProfile)
    {
        if(lineSet instanceof StsPreStackLineSet2d)
            currentLine2d = (StsPreStackLine2d)lineSet.getDataLine(velocityProfile.row, velocityProfile.col);
        else
            currentLine2d = null;
    }

    private void initializePicks()
	{
		vertex = null;
        velocityProfile.changeType = CHANGE_NONE;
//        vertexChanged = false;
    }

	public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel)
	{
        StsView currentView = glPanel.getView();

        StsVelocityProfile currentVelocityProfile = velocityVolume.getCurrentVelocityProfile(model.win3d);
        if(currentVelocityProfile == null) return false;

        if(currentVelocityProfile != velocityProfile)
        {
            checkSetCurrentLine2d(currentVelocityProfile);
            if(velocityProfile != null && !velocityProfile.sameRowCol(currentVelocityProfile))
            {
                velocityProfile.checkCurrentVelocityProfile();
                if(debug) System.out.println("new profile at row " + currentVelocityProfile.row + " col " + currentVelocityProfile.col);
                velocityProfile = currentVelocityProfile;
                initializePicks();
            }
            else
            {
                if(velocityProfile != null) System.out.println("velocityProfile && currentVelocityProfile are same row col but different objects.");
                velocityProfile = currentVelocityProfile;
            }
        }

        if(debug) System.out.println("current profile at row " + velocityProfile.row + " col " + velocityProfile.col);

		if(currentView instanceof StsViewGather)
			return handleAction(mouse, glPanel, (StsViewGather)currentView);
		else if(currentView instanceof StsView3d)
			return handleAction(mouse, glPanel, (StsView3d)currentView);
		else if(currentView instanceof StsViewCursor)
			return handleAction(mouse, glPanel, (StsViewCursor)currentView);
        else if(currentView instanceof StsViewPreStack)
            return handleAction(mouse, glPanel, (StsViewPreStack)currentView);
        return true;
	}

	private boolean handleAction(StsMouse mouse, StsGLPanel glPanel, StsView3d view3d)
	{
		if(lineSet instanceof StsPreStackLineSet3d)
            view3d.moveCursor3d(mouse, (StsGLPanel3d)glPanel);
        else
            ((StsPreStackLineSet2d)lineSet).performMouseAction(mouse, glPanel);
        return true;
	}

	private boolean handleAction(StsMouse mouse, StsGLPanel glPanel, StsViewCursor viewCursor)
	{
		return viewCursor.moveCursor3d(mouse, (StsGLPanel3d)glPanel);
	}

	private boolean handleAction(StsMouse mouse, StsGLPanel glPanel, StsViewPreStack viewSeismic)
	{
        viewSeismic.setVelocityProfile(velocityProfile);
        pick = viewSeismic.computePickPoint(mouse);
		if(mouse.getCurrentButton() == StsMouse.LEFT)
		{
            if(debug) printDebug(viewSeismic, mouse);
            int buttonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);
            if(buttonState == StsMouse.PRESSED)
			{

                StsMethod pickMethod = new StsMethod(viewSeismic, "pickOnSeismicView", new Object[] { glPanel, velocityProfile } );
				StsJOGLPick.pick3d(glPanel, pickMethod, StsMethodPick.PICKSIZE_EXTRA_LARGE, StsMethodPick.PICK_ALL);
				if(!processSeismicPicks(viewSeismic))return true;
				viewSeismic.viewObjectRepaint(viewSeismic, velocityProfile);  // redraw just this view as vertex might be highlighted
			}
			// Rubber band from last or between two on either side
			else if(buttonState == StsMouse.DRAGGED && vertex != null)
			{
				moveVertex(pick);
                model.viewObjectRepaint(viewSeismic, velocityProfile); // redraws profile as we drag it
			}
			else if(buttonState == StsMouse.RELEASED)
			{
                if(vertex != null && velocityProfile.changeType != CHANGE_NONE)
				{
                    /*         //this is no longer valid - we're letting the user pick outside viewport but snapping it to t=0 or t=traceLength SWC 8/5/09
					if(!viewSeismic.isInsideInsetViewPort(mouse.mousePoint.x, mouse.mousePoint.y))
						viewSeismic.restoreVertexPick(velocityProfile);
                     */
                    if(!viewSeismic.isInsideInsetViewPort(mouse.mousePoint.x, mouse.mousePoint.y))    //wrw reenable
                        viewSeismic.restoreVertexPick(velocityProfile);
                    checkAddProfileAndPick();

                    model.viewObjectChanged(viewSeismic, velocityProfile);
                    model.viewObjectChanged(viewSeismic, currentLine2d);
                    model.viewObjectRepaint(viewSeismic, velocityVolume); // changes 3d view of 2d data
                    model.viewObjectRepaint(viewSeismic, lineSet); // changes 3d view of 3d data (lineSet is cursor3d displayable object)
 					model.viewObjectRepaint(viewSeismic, velocityProfile); // redraws 2d line in 3d view
				}
                velocityProfile.clearInitialProfilePoints();
                indexPicked = StsGather.NONE_INDEX;
                typePicked = TYPE_NONE;
            }
		}
        return true;
	}

    private void printDebug(StsViewPreStack viewSeismic, StsMouse mouse)
    {
        System.out.println("handleAction(" + viewSeismic.getClass().getName() + ") called for " + mouse.toString());
    }

	private void moveVertex(StsPoint pick)
	{
        if(vertex.sameAs(pick)) return;
        setChangeType();

        checkVertexLimits(pick);
		if(debug) System.out.println("vertex " + vertex.toString() + " moved to: " + pick.toString());
		vertex.setValues(pick);
        velocityProfile.moveVertex( indexPicked, vertex);
	}

    /** If profile is not persistent, add it (which will include adding any existing points).
     *  If persistent, we just need to add a point.
     *  Update interpolation with new profile and point.
     * @return
     */
    private void checkAddProfileAndPick()
    {
        // if this profile is currently interpolated, make it an actual profile and add to model
        // this will also add the pick point that was added creating this profile
        if(velocityProfile.isInterpolated())
        {
            velocityVolume.addProfile(velocityProfile);
            velocityVolume.checkRunInterpolation();
        }
        // otherwise it's already an actual profile, and we just add the pick point
        else
        {
            velocityProfile.addPickTransactionCmd(typePicked, indexPicked);
            if(indexPicked >= 0)
                velocityVolume.updateVelocityProfile(velocityProfile, indexPicked);
        }
        // indexPicked = StsGather.NONE_INDEX;
    }

    private boolean processSeismicPicks(StsViewPreStack viewSeismic)
	{
		if(StsJOGLPick.hits == 0) // we have added a point
		{
            if(pick == null) return false;
            velocityProfile.clearInitialProfilePoints();
            velocityProfile.getSetVvsInitialProfilePoints();

            if(!addSemblancePoint(pick)) return false;
//            velocityProfile.setInitialProfilePoints(indexPicked, pick);
            velocityProfile.changeType = CHANGE_POINT;

//            vertexChanged = true;
            return true;
        }
		StsPickItem pickItem = StsJOGLPick.pickItems[0];
		typePicked = pickItem.names[0];
		indexPicked = pickItem.names[1];
		if(typePicked == TYPE_GATHER_VERTEX) // we have selected a point
		{
			if(debug) System.out.println("vertex " + indexPicked + " picked: " + pick.toString());
            velocityProfile.setInitialProfilePoints(indexPicked, pick);
            velocityProfile.getSetVvsInitialProfilePoints();
            vertex = velocityProfile.getProfilePoint(indexPicked);
//            velocityProfile.setPickPoint(indexPicked);
            viewSeismic.setVertexPick(indexPicked, velocityProfile);
            return true;
		}
		else if(typePicked == TYPE_GATHER_EDGE && indexPicked >= 0) // we've added a point to an edge
		{
            viewSeismic.setVertexPick(-1, velocityProfile);
			if(debug) System.out.println("edge segment " + indexPicked + " picked. Vertex created: " + pick.toString());
            if(!addSemblancePoint(pick)) return false;
//            velocityProfile.setInitialProfilePoints(indexPicked, pick);
//            velocityProfile.setPickPoint(indexPicked);
            return true;
        }
		else
		{
			StsException.systemError("StsVelocityAnalysisEdit.processSeismicPicks() failed. Found a pick but couldn't identify.");
            velocityProfile.clearInitialProfilePoints();
//            velocityProfile.clearPickPoint();
            return false;
		}
	}

    /** If this is a new profile (not peristed yet), it probably has some interpolated points which we will now clear.
     *  We then add this first pick point, but won't commit this until the mouse is released.
     *  If persisted already, it has at least one point.
     * @param pick point picked
     * @return unconditionally true (add catch and return false if this is a problem)
     */
    private boolean addSemblancePoint(StsPoint pick)
    {
        if(pick == null) return false;
        if(velocityProfile.isInterpolated())
            velocityProfile.clearProfilePoints();
        vertex = pick;
        indexPicked = velocityProfile.addProfilePoint(pick);
        typePicked = velocityProfile.typePicked;
        return true;
    }

	private boolean processGatherPicks(StsViewGather viewGather)
	{
		if(StsJOGLPick.hits == 0)
        {
            indexPicked = StsGather.NONE_INDEX;
            typePicked = TYPE_NONE;
            return true;
        }
		StsPickItem pickItem = StsJOGLPick.pickItems[0];
		typePicked = pickItem.names[0];
		indexPicked = pickItem.names[1];

        if(indexPicked == StsGather.STRETCH_MUTE_INDEX)
		{
			if(velocityProfile.getProfilePoints() == null)
			{
				new StsMessage(viewGather.getWindow(), StsMessage.WARNING, "Can't adjust stretch mute until a velocity profile has been picked on semblanceBytes.", true);
				return true;
			}
			vertex = pick;
		}
        else // top or bottom mute pick or NMO curve pick
        {
            vertex = velocityProfile.getProfilePoint(indexPicked);
            velocityProfile.setInitialProfilePoints(indexPicked, pick);
        }

        if(debug)
		{
			if (typePicked == TYPE_GATHER_VERTEX)
				System.out.println("gather vertex " + indexPicked + " picked: " + vertex.toString());
			else if(indexPicked >= 0)
				System.out.println("gather edge " + indexPicked + " picked: " + pick.toString());
			else if(indexPicked >= StsGather.TOP_MUTE_INDEX)
				System.out.println("gather edge top mute picked: " + pick.toString());
			else if(indexPicked >= StsGather.BOT_MUTE_INDEX)
				System.out.println("gather edge bottom mute picked: " + pick.toString());
			else if(indexPicked >= StsGather.STRETCH_MUTE_INDEX)
				System.out.println("gather edge stretch mute picked: " + pick.toString());
		}
		return true;
	}

    private void setChangeType()
    {
        if(indexPicked < 0)
            velocityProfile.changeType = StsVelocityProfile.CHANGE_MUTE;
        else
            velocityProfile.changeType = StsVelocityProfile.CHANGE_POINT;
    }

    private void moveVertexOnGather(StsPoint pick, StsSuperGather gather)
	{
        setChangeType();
        if(indexPicked == StsGather.STRETCH_MUTE_INDEX)
        {
            gather.adjustStretchMute(vertex, pick);
        }
        else if(this.typePicked == TYPE_GATHER_EDGE)
		{
			velocityProfile.adjustVelocityOnGather(vertex, indexPicked, pick, gather);
			checkVertexVelocityLimit(vertex);
            velocityProfile.moveVertex( indexPicked, vertex);
        }
		else if(typePicked == TYPE_GATHER_VERTEX)
		{
			velocityProfile.adjustTzeroOnGather(vertex, pick);
			checkVertexTimeLimit(vertex);
            velocityProfile.moveVertex( indexPicked, vertex);
        }
    }

    private void checkVertexLimits(StsPoint pick)
	{
        if(indexPicked < 0)
        {
            return;
        }

        if(indexPicked >= 1)
		{
			StsPoint pointAbove = velocityProfile.getProfilePoint(indexPicked - 1);
			float v0 = pointAbove.v[0];
			float t0 = pointAbove.v[1];
			float t1 = pick.v[1];
			if(t1 <= t0)
			{
				t1 = t0 + 1;
				pick.v[1] = t1;
				StsToolkit.beep();
			}
			float v1 = pick.v[0];
			if(v1 * v1 * t1 < v0 * v0 * t0)
			{
				pick.v[0] = v0 + 0.01f;
				StsToolkit.beep();
			}
		}
		if(indexPicked < velocityProfile.getNProfilePoints() - 1)
		{
			StsPoint pointBelow = velocityProfile.getProfilePoint(indexPicked + 1);
			float v1 = pointBelow.v[0];
			float t1 = pointBelow.v[1];
			float t0 = pick.v[1];
			if(t0 >= t1)
			{
				t0 = t1 - 1;
				pick.v[1] = t0;
				StsToolkit.beep();
			}
			float v0 = pick.v[0];
			if(v0 * v0 * t0 > v1 * v1 * t1)
			{
				pick.v[0] = v1 - 0.01f;
				StsToolkit.beep();
			}
		}
		if (pick.v[1] < lineSet.zMin || pick.v[1] > lineSet.zMax)
		{
		    pick.v[1] = Math.max(lineSet.zMin, pick.v[1]);
		    pick.v[1] = Math.min(lineSet.zMax, pick.v[1]);
		    StsToolkit.beep();
		}
	}

	private void checkVertexTimeLimit(StsPoint pick)
	{
        if(indexPicked < 0)
        {
            return;
        }
        if(indexPicked >= 1)
		{
			StsPoint pointAbove = velocityProfile.getProfilePoint(indexPicked - 1);
			float t0 = pointAbove.v[1];
			float t1 = pick.v[1];
			if(t1 <= t0)
			{
				t1 = t0 + 1;
				pick.v[1] = t1;
				StsToolkit.beep();
			}
		}
		if(indexPicked < velocityProfile.getNProfilePoints() - 1)
		{
			StsPoint pointBelow = velocityProfile.getProfilePoint(indexPicked + 1);
			float t1 = pointBelow.v[1];
			float t0 = pick.v[1];
			if(t0 >= t1)
			{
				t0 = t1 - 1;
				pick.v[1] = t0;
				StsToolkit.beep();
			}
		}
	}

	private void checkVertexVelocityLimit(StsPoint pick)
	{
        if(indexPicked < 0) return;

        if(indexPicked >= 1)
		{
			StsPoint pointAbove = velocityProfile.getProfilePoint(indexPicked - 1);
			float v0 = pointAbove.v[0];
			float t0 = pointAbove.v[1];
			float t1 = pick.v[1];
			float v1 = pick.v[0];
			if(v1 * v1 * t1 < v0 * v0 * t0)
			{
				pick.v[0] = v0 + 0.01f;
				StsToolkit.beep();
			}
		}
		if(indexPicked < velocityProfile.getNProfilePoints() - 1)
		{
			StsPoint pointBelow = velocityProfile.getProfilePoint(indexPicked + 1);
			float v1 = pointBelow.v[0];
			float t1 = pointBelow.v[1];
			float t0 = pick.v[1];
			float v0 = pick.v[0];
			if(v0 * v0 * t0 > v1 * v1 * t1)
			{
				pick.v[0] = v1 - 0.01f;
				StsToolkit.beep();
			}
		}
	}

	//Replacing key action with dialog on vertex or line selection.  TJL 7/7/07
    public boolean keyReleased(KeyEvent e, StsMouse mouse, StsGLPanel glPanel)
    {
        if(mouse.isButtonDown()) return false;

        if(e.getKeyCode() == KeyEvent.VK_D)
        {
            if (velocityProfile == null) return false;
            StsVelocityCursor cursor = StsVelocityCursor.getStsVelocityCursor();
            if (cursor == null) return false;
            double time = cursor.getTime();
            StsPoint point = velocityProfile.getClosestProfilePoint(time);
            if (point == null) return false;
            velocityProfile.deleteProfilePoint(point);
            model.viewObjectRepaint(this, velocityProfile);
            return true;
        }
        else
            return false;
    }

	private boolean handleAction(StsMouse mouse, StsGLPanel glPanel, StsViewGather viewGather)
	{
		int buttonState;

		pick = viewGather.computePickPoint(mouse);

		if(mouse.getCurrentButton() == StsMouse.LEFT)
		{
			buttonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);
			if(debug) mouse.printState();
			if(buttonState == StsMouse.PRESSED)
			{
				StsMethod pickMethod = new StsMethod(velocityProfile, "pickOnGather", new Object[] { glPanel, (StsViewGather)viewGather} );
				//StsMethodPick picker = new StsMethodPick(glPanel, pickMethod, StsMethodPick.PICKSIZE_MEDIUM, StsMethodPick.PICK_CLOSEST);
				StsJOGLPick picker = new StsJOGLPick(glPanel, pickMethod, StsMethodPick.PICKSIZE_MEDIUM, StsMethodPick.PICK_ALL);
				try
				{
					//glPanel3d.glc.pick(picker);
					picker.methodPick3d();
					if(!processGatherPicks(viewGather)) return true;
					if(debug) System.out.println("Velocity curve " + nPickedCurve + " picked." + " coordinates " + pick.v[0] + " " + pick.v[1]);
				}
				catch(Exception e)
				{
					StsException.outputException("StsPolygonAction.handlePolygonAction() failed.", e, StsException.WARNING);
				}
				viewGather.checkSetOrder(StsSuperGather.ORDER_2ND);
			}
			// Rubber band from last or between two on either side
			else if (buttonState == StsMouse.DRAGGED && vertex != null)
			{
				moveVertexOnGather(pick, viewGather.superGather);
//				viewGather.superGather.setUseSecondOrder(true);
                viewGather.viewObjectRepaint(this, velocityProfile); // redraw only this view while we drag
//                model.viewObjectRepaint(this, velocityProfile);
            }
			else if (buttonState == StsMouse.RELEASED) // add a point unless near firstPoint && nPicks > 3
			{
                if(indexPicked == StsGather.NONE_INDEX) return true;
                viewGather.superGather.setUseSecondOrder(false);
                byte order = lineSet.getSemblanceComputeProperties().order;
				viewGather.checkSetOrder(order);
                if(indexPicked >= 0)
                    checkAddProfileAndPick();
//                if(indexPicked >= 0)
//                    velocityVolume.updateVelocityProfile(velocityProfile, indexPicked);
                if(indexPicked == StsGather.STRETCH_MUTE_INDEX)
                    stretchMuteChanged();
                model.viewObjectChanged(this, velocityProfile);
                model.viewObjectChanged(this, currentLine2d);
//                velocityProfile.changeType = CHANGE_NONE;
                model.viewObjectRepaint(this, velocityVolume); // changes 3d view of 2d velocity data
                model.viewObjectRepaint(this, lineSet); // changes 3d view of 2d stacked data
                model.viewObjectRepaint(this, velocityProfile); // redraws 2d line in 3d view
                // velocityProfile.clearInitialProfilePoints();
                indexPicked = StsGather.NONE_INDEX;
                typePicked = TYPE_NONE;
            }
        }
		else if(mouse.getCurrentButton() == StsMouse.VIEW)
		{
			buttonState = mouse.getButtonStateCheckClear(StsMouse.VIEW);

			// Delete Vertex
			if (buttonState == StsMouse.RELEASED)
			{
                viewGather.viewObjectRepaint(this, velocityProfile);
			}
		}
//		model.win3dDisplayAll();
		return true;
	}

    private boolean stretchMuteChanged()
    {
        StsWiggleDisplayProperties wiggleProperties = lineSet.getWiggleDisplayProperties();
        if(!wiggleProperties.stretchMuteChanged) return false;
        return model.viewObjectChanged(this, wiggleProperties);
    }
/*
    private boolean checkCommitVelocityProfile()
    {
        if(!velocityProfile.addProfileOk()) return false;
        endAndStartTransaction();
        return true;
    }
*/
/*
    private void endAndStartTransaction()
	{
		velocityVolume.commitCurrentVelocityProfile();
		actionManager.commit();
		actionManager.getCreateCurrentTransaction("StsVelocityAnalysisEdit");
	}
*/
	public boolean end()
	{
		if(velocityProfile != null) velocityProfile.checkCurrentVelocityProfile();
		lineSet.setMode(StsPreStackLineSet.DISPLAY_MODE);
		if (statusArea != null)
		{
			statusArea.textOnly(); //sometimes this is null??? SWC 5/4/09
		}
        StsToggleButton toggleButton = (StsToggleButton)model.win3d.getToolbarComponentNamed(StsVelocityAnalysisToolbar.NAME,  StsVelocityAnalysisToolbar.EDIT_VELOCITY);
        toggleButton.setSelected(false);
//		model.win3dDisplayAll();  // redraw in case display&edit modes are different
		return true;
	}

    public boolean canBeInterrupted()
    {
        if(velocityVolume != null)
            velocityVolume.checkCurrentVelocityProfile(model.win3d);
        model.commit();
        return true;
    }
}
