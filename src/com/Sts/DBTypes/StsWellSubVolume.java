package com.Sts.DBTypes;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */

import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;

public class StsWellSubVolume extends StsObject
{
    StsGridPoint svCenter = null;
    StsSeismicVolume volume;
    StsWell well;
    float top = 10.0f;     // offset from center along wellbore
    float bottom = 10.0f;  // offset from center along wellbore
    float radius = 10.0f;  // radius of cylinder

    public StsWellSubVolume()
    {
    }

    public StsWellSubVolume(boolean persistent)
    {
        super(persistent);
    }

    public StsWellSubVolume(StsCursorPoint cursorPoint, StsWell well, float topOffset, float btmOffset, float radius)
    {
        this.volume = null;
        this.well = well;
        float[] xyz = cursorPoint.point.v;
        svCenter = new StsGridPoint(cursorPoint.point, volume);
        setTop(topOffset);
        setBtm(btmOffset);
        setRadius(radius);
    }

    public StsWellSubVolume(StsCursorPoint cursorPoint, StsSeismicVolume volume, StsWell well, float topOffset, float btmOffset, float radius)
    {
        this.volume = volume;
        this.well = well;
        float[] xyz = cursorPoint.point.v;
        svCenter = new StsGridPoint(cursorPoint.point, volume);
        setTop(topOffset);
        setBtm(btmOffset);
        setRadius(radius);
    }

    public StsSeismicVolume getVolume()
    {
        return volume;
    }

    public void display(StsGLPanel glPanel, StsColor stsColor, byte action, boolean editing)
    {
        GL gl = glPanel.getGL();
        displayBoundingCylinder(gl, stsColor, StsGraphicParameters.gridLineWidth);

        if (svCenter == null)
        {
            return;
        }

        int boxSize = 4;
        if (editing && action == StsWellSetSubVolume.ACTION_MOVE_CYL)
            boxSize = 8;

        StsGLDraw.drawPoint(svCenter.getXYZorT(), StsColor.BLACK, glPanel, boxSize + 4, 1.0);
        StsGLDraw.drawPoint(svCenter.getXYZorT(), stsColor, glPanel, boxSize, 2.0);
        displayBoundingCylinder(gl, stsColor, StsGraphicParameters.gridLineWidth);
        if (editing && action == StsBoxSetSubVolume.ACTION_MOVE_POINT)
            drawFacePoints(gl, false);
    }

    public void pickCenterPoint(GL gl, boolean editing, int boxID)
    {
        if (svCenter == null)
        {
            return;
        }

        int boxSize = 4;
        if (editing) boxSize = 8;

        gl.glInitNames();
        gl.glPushName(boxID);
        gl.glPushName( -1);
        StsGLDraw.drawPoint(svCenter.getXYZorT(), gl, boxSize + 4);
        gl.glPopName();
        gl.glPopName();
    }

    public void pickFacePoint(GL gl, int cylID)
    {
        gl.glInitNames();
        gl.glPushName(cylID);
        drawFacePoints(gl, true);
    }

    public void moveCylinder(StsMouse mouse, StsGLPanel3d glPanel3d)
    {
        StsPoint newCenterPoint;

        StsView currentView = glPanel3d.getView();
        if (currentView instanceof StsView3d)
        {
            int planeDir;

            double[] lineVector = glPanel3d.getViewVectorAtMouse(mouse);
            if (Math.abs(lineVector[0]) > Math.abs(lineVector[1]))
            {
                planeDir = 0;
                if (Math.abs(lineVector[2]) >= Math.abs(lineVector[0]))
                    planeDir = 2;
            }
            else
            {
                planeDir = 1;
                if (Math.abs(lineVector[2]) >= Math.abs(lineVector[1]))
                    planeDir = 2;
            }
            StsCursor3d cursor3d = glPanel3d.getCursor3d();
            newCenterPoint = cursor3d.getPointInPlaneAtMouse(glPanel3d, planeDir, svCenter.getPoint(), mouse);
            if (newCenterPoint == null)
            {
                return;
            }
        }
        else if (currentView instanceof StsViewCursor)
        {
            StsCursorPoint cubePoint = ( (StsViewCursor) currentView).getCursorPoint(mouse);
            if (cubePoint == null)
            {
                return;
            }
            newCenterPoint = cubePoint.point;
        }
        else
        {
            return;
        }

        StsPoint dCenter = StsPoint.subPointsStatic(newCenterPoint, svCenter.getPoint());
        adjustXYZPosition(dCenter);
        svCenter = new StsGridPoint(newCenterPoint, volume);
    }

    public void movePoint(StsMouse mouse, StsGLPanel3d glPanel3d, int pickedPointIndex)
    {
        float[] faceCenter;
        StsPoint newFaceCenterPoint;

        StsView currentView = glPanel3d.getView();
        if (! (currentView instanceof StsView3d))
        {
            return;
        }

        // X and Y face points will be moved in the Z-plane thru the face point
        // Z face points will be moved in the X or Y planes depending on orientation
        int planeDir = 2;
        if (pickedPointIndex > 3)
        {
            double[] lineVector = glPanel3d.getViewVectorAtMouse(mouse);
            if (Math.abs(lineVector[0]) > Math.abs(lineVector[1]))
            {
                planeDir = 0;
            }
            else
            {
                planeDir = 1;
            }
        }
        int moveDir = pickedPointIndex / 2;
        faceCenter = getFaceCenter(pickedPointIndex);
        if(faceCenter == null) return;
        StsCursor3d cursor3d = glPanel3d.getCursor3d();
        newFaceCenterPoint = cursor3d.getPointInPlaneAtMouse(glPanel3d, planeDir, faceCenter, mouse);
        if (newFaceCenterPoint == null) return;
        float[] dCenter = StsMath.subtract(newFaceCenterPoint.v, faceCenter);
        adjustRange(pickedPointIndex, dCenter[moveDir]);
        resetCenter();
    }

    private void resetCenter()
    {
        // Add center reset code.....deviated wells will be interesting. Reset along wellbore path.
    }

    private void adjustRange(int pickPt, float delta)
    {
        // Add adjust cylinder code....
    }

    private void adjustXYZPosition(StsPoint position)
    {
        // Add adjust center of cylinder code along wellbore path....
    }

    private float[] getFaceCenter(int pickedPointIndex)
    {
        // Add code to get face center for top or bottom...Face is always perpendicular to wellbore
        return null;
    }

    public void displayBoundingCylinder(GL gl, StsColor stsColor, float lineWidth)
    {
        // Add code to display a wireframe of the bounding cylinder...
    }

    public void drawFacePoints(GL gl, boolean isPicking)
    {
        // Add code to display the points on the top and bottom faces
    }

    public void setTop(float top) { this.top = top; }   // Top offset from center along the wellbore
    public void setBtm(float btm) { this.bottom = btm; }  // Bottom offset from center along wellbore
    public void setRadius(float radius) { this.radius = radius; }

    public float getTop() { return top; }
    public float getBtm() { return bottom; }
    public float getRadius() { return radius; }

}
