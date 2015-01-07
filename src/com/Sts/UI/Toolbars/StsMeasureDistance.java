package com.Sts.UI.Toolbars;

import com.Sts.Actions.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Reflect.StsMethod;
import com.Sts.Types.*;
import com.Sts.UI.StsMessage;
import com.Sts.Utilities.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Mar 8, 2008
 * Time: 2:59:08 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsMeasureDistance extends StsAction
{
    StsView currentView;
    PickedObject firstPick;
    PickedObject lastPick;
    int rulerType = StsMouseActionToolbar.NONE;
    StsLineMeasure connectingLine = null;
	double[] screenPoint = new double[3];

	static final int pixelRadius = 10;

    public StsMeasureDistance(StsActionManager actionManager)
    {
        super(actionManager, true);

        if(actionManager.getCurrentAction() instanceof StsMeasureDistance)
        {
            actionManager.endCurrentAction();
        }
        if(actionManager.getModel().getProject().getZDomain() != StsProject.TD_DEPTH)
        {
            new StsMessage(actionManager.getModel().win3d, StsMessage.WARNING, "Can only use measuring device in Depth.\n Switch domains and try again.");
            end();
        }
    }

    public boolean start()
    {
        connectingLine = new StsLineMeasure();
        model.addDisplayableInstance(connectingLine);
        rulerType = ((StsMouseActionToolbar)getModel().win3d.getToolbarNamed(StsMouseActionToolbar.NAME)).getRulerType();
        return true;
    }

    public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel)
    {
        if(rulerType == StsMouseActionToolbar.NONE)
            return false;

        StsView view = ((StsGLPanel3d)glPanel).getView();
        if(view instanceof StsView3d)
        {
            if(currentView != view)
            {
                clearObjects();
                currentView = view;
            }

            if(rulerType == StsMouseActionToolbar.BETWEEN_OBJECTS)
                return measureBetweenObjects(view, mouse, glPanel);
            else
                return  measureBetweenPoints(view, mouse, glPanel);
        }
        else
        {
            StsMessageFiles.infoMessage("Ruler operation only works in 3D view.");
            end();
            return false;
        }
    }

    public boolean drawInitialRuler(StsMouse mouse, StsGLPanel3d glPanel)
    {
		double[] centerViewPoint = ((StsView3d)currentView).getCenterViewPoint().getPointXYZDoubles();
        double[] centerScreenPoint = glPanel.getScreenCoordinates(centerViewPoint);
		screenPoint[2] = centerScreenPoint[2];

        StsMousePoint mousePointGL = glPanel.getMousePointGL();
        screenPoint[0] = (double)mousePointGL.x-100;
        screenPoint[1] = (double)mousePointGL.y;
        connectingLine.setFirstPoint(glPanel.getWorldCoordinates(screenPoint));
        screenPoint[0] = (double)mousePointGL.x+100;
        connectingLine.setLastPoint(glPanel.getWorldCoordinates(screenPoint));

		screenPoint[0] = centerScreenPoint[0] + pixelRadius;
		screenPoint[1] = centerScreenPoint[1];
		double[] edgePoint = glPanel.getWorldCoordinates(screenPoint);
		connectingLine.radius = (float)StsMath.distance(centerViewPoint, edgePoint, 3);
        connectingLine.setPlacement(connectingLine.END);
        connectingLine.setAddEnds(true);
        model.viewObjectRepaint(this);
        return true;
    }

    public boolean measureBetweenPoints(StsView view, StsMouse mouse, StsGLPanel glPanel)
    {
       int mouseButtonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);
       StsGLPanel3d glPanel3d = (StsGLPanel3d)glPanel;

       if(mouseButtonState == StsMouse.PRESSED)
        {
            PickedObject pick =  getPickedEnd(glPanel3d, mouse);
            if(pick == null)
            {
                drawInitialRuler(mouse, glPanel3d);
                firstPick = null;
                return true;
            }
            firstPick = pick;
        }
        else if((mouseButtonState == StsMouse.DRAGGED) && (firstPick != null))
        {
            double[] screenPt;
            if(firstPick.index == connectingLine.START)
            {
                screenPt = glPanel3d.getScreenCoordinates(connectingLine.getFirstPoint());
                StsMousePoint mousePointGL = glPanel.getMousePointGL();
                screenPt[0] = (double)mousePointGL.x;
                screenPt[1] = (double)mousePointGL.y;
                connectingLine.setFirstPoint(glPanel.getWorldCoordinates(screenPt));
            }
            else
            {
                screenPt = glPanel3d.getScreenCoordinates(connectingLine.getLastPoint());
                StsMousePoint mousePointGL = glPanel.getMousePointGL();
                screenPt[0] = (double)mousePointGL.x;
                screenPt[1] = (double)mousePointGL.y;
                connectingLine.setLastPoint(glPanel.getWorldCoordinates(screenPt));
            }
            annotateLine();            
            model.viewObjectRepaint(this);
        }
        model.viewObjectRepaint(this);
        return true;
    }

    public boolean measureBetweenObjects(StsView view, StsMouse mouse, StsGLPanel glPanel)
    {
       int mouseButtonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);
       StsGLPanel3d glPanel3d = (StsGLPanel3d)glPanel;

       if(mouseButtonState == StsMouse.PRESSED && firstPick == null)
        {
            PickedObject pick = getPickedObject(glPanel3d, mouse, "First object: ");
            if(pick == null) return true;
            if(pick.sameAs(firstPick)) return true;
            firstPick = pick;
            StsMessageFiles.infoMessage("Picked object: " + firstPick.description);
            connectingLine.setFirstPoint(firstPick.point);
        }
        else if(mouseButtonState == StsMouse.DRAGGED && firstPick != null)
        {
            PickedObject pick = getPickedObject(glPanel3d, mouse, "Last object: ");
            if(pick == null)
            {
                return true;
            }
            if(pick.sameAs(lastPick) || pick.sameAs(firstPick))
            {
                return true;
            }
            lastPick = pick;
            connectingLine.setLastPoint(lastPick.point);

            annotateLine();
            model.viewObjectRepaint(this, firstPick.object);
        }
        else if(mouseButtonState == StsMouse.RELEASED)
        {
            clearObjects();
        }
        return true;
    }

    private void annotateLine()
    {
        float distance = connectingLine.getDistance();
	    float projectAngle = model.getProject().getAngle();
        double ccwFromEast = connectingLine.getAngle() + projectAngle;
		double azimuth = StsMath.getAzimuthFromCCWEast(ccwFromEast); // cwFromNorth
        double dip = Math.toDegrees(connectingLine.getDip());

        connectingLine.setShowLabel(true);
        connectingLine.setPlacement(StsLineMeasure.END);
        if(rulerType == StsMouseActionToolbar.FREE_SPACE)
            connectingLine.setAddEnds(true);
        else
            connectingLine.setAddEnds(false);

        StsDecimalFormat floatFormat = new StsDecimalFormat(6);

        connectingLine.setLabel(new String("   " + floatFormat.formatValue(distance) + " " + model.getProject().getXyUnitString() +
            ", Azimuth=" + floatFormat.formatValue(azimuth) + ", Dip=" + floatFormat.formatValue(dip)));

        StsMessageFiles.infoMessage("  " + floatFormat.formatValue(distance) + " " + model.getProject().getXyUnitString() + " between selected objects is." +
            		" Vertical=" + floatFormat.formatValue(connectingLine.getVerticalDistance()) +
            		" Horizontal=" + floatFormat.formatValue(connectingLine.getHorizontalDistance()) +
                    " Azimuth=" + floatFormat.formatValue(azimuth) + " Dip=" + floatFormat.formatValue(dip));
    }

    private PickedObject getPickedEnd(StsGLPanel3d glPanel, StsMouse mouse)
    {
            StsMethod pickMethod = new StsMethod(connectingLine, "pickPoint", new Object[] { glPanel} );
			StsJOGLPick picker = new StsJOGLPick(glPanel, pickMethod, StsMethodPick.PICKSIZE_MEDIUM, StsMethodPick.PICK_FIRST);
			if(picker.methodPick3d())
            {
                if(StsJOGLPick.pickItems[0].names[0] == connectingLine.START)
                    return new PickedObject(null, connectingLine.getPoint(StsJOGLPick.pickItems[0].names[0]), StsJOGLPick.pickItems[0].names[0], "First Point");
                else if(StsJOGLPick.pickItems[0].names[0] == connectingLine.END)
                    return new PickedObject(null, connectingLine.getPoint(StsJOGLPick.pickItems[0].names[0]), StsJOGLPick.pickItems[0].names[0], "Last Point");
            }
            return null;
    }

     private PickedObject getPickedObject(StsGLPanel3d glPanel3d, StsMouse mouse, String message)
     {
         int sensorIdx = 0;
         StsDynamicSensor sensor = (StsDynamicSensor)StsJOGLPick.pickVisibleClass3d(glPanel3d, StsDynamicSensor.class, StsJOGLPick.PICKSIZE_MEDIUM, StsJOGLPick.PICK_FIRST);
         if(sensor != null)
         {
            StsPickItem items = StsJOGLPick.pickItems[0];
            if(items.names.length < 2) return new PickedObject(sensor, null, 0);
            sensorIdx = items.names[1];
            StsPoint point = new StsPoint(sensor.getXYZMT(sensorIdx));
            return new PickedObject(sensor, point, sensorIdx);
         }
         StsWell well = (StsWell) StsJOGLPick.pickVisibleClass3d(glPanel3d, StsWell.class, StsJOGLPick.PICKSIZE_MEDIUM, StsJOGLPick.PICK_CLOSEST);
         if (well != null)
         {
            StsPickItem item = StsJOGLPick.pickItems[0];
		    int nSegment = item.names[1];
            StsPoint point = well.getPointOnLineNearestMouse(nSegment, mouse, glPanel3d);
            return new PickedObject(well, point, nSegment);
         }
         return null;
     }

    public void clearObjects()
    {
        firstPick = null;
        lastPick = null;
        connectingLine.clearPoints();
    }

    public boolean end()
    {
        model.removeDisplayableInstance(connectingLine);
        connectingLine = null;
        model.viewObjectRepaint(this);
        return true;
    }

    class PickedObject
    {
        StsObject object;
        StsPoint point;
        int index;
        String description;

        PickedObject(StsObject object, StsPoint point, int index)
        {
            this.object = object;
            this.point = point;
            this.index = index;
            if(object != null)
                description = object.toString();
        }

        PickedObject(StsObject object, StsPoint point, int index, String description)
        {
            this.object = object;
            this.point = point;
            this.index = index;
            this.description = description;
        }

        PickedObject(StsSensor sensor, StsPoint point, int index)
        {
            this.object = sensor;
            this.point = point;
            this.index = index;
            description = sensor.toString(index);
        }

        boolean sameAs(PickedObject otherObject)
        {
            if(otherObject == null) return false;
            return otherObject.object == object && otherObject.index == index;
        }
    }
}