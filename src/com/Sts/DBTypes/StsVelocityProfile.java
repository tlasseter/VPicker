package com.Sts.DBTypes;

import com.Sts.Actions.Wizards.VelocityAnalysis.*;
import com.Sts.DB.DBCommand.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.Types.PreStack.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.io.*;
import java.text.*;
import java.util.*;

/**
 * <p>Title: S2S development</p>
 * <p/>
 * <p>Description: Integrated seismic to simulation software</p>
 * <p/>
 * <p>Copyright: Copyright (c) 2004</p>
 * <p/>
 * <p>Company: S2S Systems LLC</p>
 *
 * @author not attributable
 * @version c51c
 */
public class StsVelocityProfile extends StsObject implements Cloneable, Serializable, Comparable
{
    /** row of gather this profile is on; in 3d this is the row number in the volume; in 2d this is the line number */
    public int row;
    /** col of gather htis profile is on; in 3d this is the col number in the volume; in 2d this is the gather sequence number on line */
    public int col;
    /** a 2d point defining the time and velocity for this top mute; mute curve is 2d-order NMO curve. trace points above this curve are muted */
    StsPoint topMute = null;
    /** a 2d point defining the time and velocity for this bottom mute; mute curve is 2d-order NMO curve. trace points below this curve are muted */
    StsPoint bottomMute = null;
    /** 2D point array of semblance time-velocity points; coordinates 0: velocity, 1: time. */
    private StsPoint[] profilePoints = new StsPoint[0];
    /** initial set of profile points before editing. Used in adjusting corridor. */
    transient public StsPoint[] initialProfilePoints;
    /** Profile set when vvs editing is started.  Maintained until vvs is updated and then deleted. */
    transient public StsPoint[] vvsInitialProfilePoints;
    /** indicates profile has been changed and model needs to be updated */
    transient public byte changeType = CHANGE_NONE;
    /** indicates profilePoints are interpolated values (not picked) */
    transient public boolean interpolated = true;
    /** index of current semblance point selected */
    transient public int indexPicked = StsGather.NONE_INDEX;
    /** current point which has been picked and may be in the process of being edited; this is the original value */
    transient public StsPoint initialPickedPoint;

    private StsPreStackLineSet lineSet = null;
    public int typePicked;

    static public final byte CHANGE_NONE = 0;
    static public final byte CHANGE_POINT = 1;
    static public final byte CHANGE_MUTE = 2;

    static final float largeFloat = StsParameters.largeFloat;
    static final boolean debug = false;


    static public final StsColor VERTEX_BORDER_COLOR = StsColor.DARK_BLUE;
    static public final StsColor VERTEX_COLOR = StsColor.PEACHPUFF;
    static public final StsColor LABEL_COLOR = StsColor.WHITE;

    static PointComparator pointComparator = new PointComparator();

    public StsVelocityProfile()
    {
    }

    public StsVelocityProfile(boolean persistent)
    {
        super(persistent);
    }

    public StsVelocityProfile(StsPreStackLineSet lineSet, float rowNum, float colNum)
    {
        super(false);
        this.lineSet = lineSet;
        this.row = lineSet.getRowFromRowNum(rowNum);
        this.col = lineSet.getColFromColNum(colNum);
        //       this.lineSetClass = (StsPreStackLineSet3dClass)lineSet.getCreateStsClass();
//		setDefaultMutes(lineSet);
    }

    public StsVelocityProfile(int row, int col)
    {
        super(false);
        lineSet = StsPreStackLineSetClass.currentProjectPreStackLineSet;
        this.row = row;
        this.col = col;
//		 this.lineSetClass = (StsPreStackLineSet3dClass)lineSet.getCreateStsClass();
//		 setDefaultMutes(lineSet);
    }

    public StsVelocityProfile(int row, int col, StsVelocityProfile otherProfile)
    {
        super(false);
        lineSet = StsPreStackLineSetClass.currentProjectPreStackLineSet;
        this.row = row;
        this.col = col;
        setProfilePoints(StsPoint.copy(otherProfile.profilePoints));
        topMute = new StsPoint(otherProfile.topMute);
        bottomMute = new StsPoint(otherProfile.bottomMute);
    }

    /*
        private void setDefaultMutes(StsPreStackLineSet lineSet)
        {
            StsSemblanceRangeProperties rangeProperties = lineSet.semblanceRangeProperties;
            float zMin = rangeProperties.getZMin();
            float zMax = rangeProperties.getZMax();
            topMute = new StsPoint(StsParameters.largeFloat, zMin);
            bottomMute = new StsPoint(StsParameters.largeFloat, zMax);
    //		stretchMute = lineSet.wiggleDisplayProperties.stretchMute;
        }
    */
    public boolean initialize(StsModel model)
    {
        int nPoints = getProfilePoints().length;
        if (nPoints <= 1)
        {
            System.out.println("Profile found with " + nPoints + " points at row " + row + " col " + col + ". Will delete.");
//            delete();
        }
        return true;
    }

    public boolean isInterpolated()
    {
        return getIndex() == -1;
    }

    public void checkCurrentVelocityProfile()
    {
        if (!isInterpolated() && getNProfilePoints() < 2)
        {
            boolean keep = StsYesNoDialog.questionValue(currentModel.win3d, "Previous profile has " + getNProfilePoints() + " points.\n" +
                "Do you wish to keep it?  If not, it will be deleted");
            if (!keep) delete();
        }
    }

    public void initializeScratchProfile(StsPreStackLineSet lineSet, int row, int col)
    {
        this.row = row;
        this.col = col;
    }

    public StsPoint getProfilePoint(int index)
    {
        if (index == StsGather.TOP_MUTE_INDEX) return getTopMute();
        else if (index == StsGather.BOT_MUTE_INDEX) return getBottomMute();
        else
        {
            StsPoint[] profilePoints = getProfilePoints();
            return profilePoints[index];
        }
    }

    public StsPoint getProfilePointCopy(int index)
    {
        if (index == StsGather.TOP_MUTE_INDEX) return getTopMute().copy();
        else if (index == StsGather.BOT_MUTE_INDEX) return getBottomMute().copy();
        else
        {
            StsPoint[] profilePoints = getProfilePoints();
            return profilePoints[index].copy();
        }
    }

    static public boolean isMuteIndex(int pointIndex)
    {
        return pointIndex == StsGather.TOP_MUTE_INDEX || pointIndex == StsGather.BOT_MUTE_INDEX;
    }

    public void setProfilePoint(int index, StsPoint point)
    {
        if (lineSet != null && lineSet.getDatumShift() != 0.0)
        {
            StsPoint pointCopy = new StsPoint(point);
            pointCopy.shiftY(-(float) lineSet.getDatumShift());
            getProfilePoints()[index] = pointCopy;
        }
        else
        {
            if (index >= getProfilePoints().length || index < 0)
                System.out.println("StsVelocityProfile: trying to set profile index: " + index + " when length is: " + getProfilePoints().length);
            else
                getProfilePoints()[index] = point;
        }
    }

    /**
     * When a point on profile has been selected and will be moved, save its original position;
     * clear when finished.
     */
/*
    public void setPickPoint(int index)
    {
        System.out.println("StsVelocityProfile.setPickPoint: " + index);
        indexPicked = index;
        initialPickedPoint = profilePoints[index].copy();
    }
*/
/*
    public void clearPickPoint()
     {
         System.out.println("StsVelocityProfile.clearPickPoint: " + indexPicked);
         indexPicked = StsGather.NONE_INDEX;
         initialPickedPoint = null;
     }
 */
    public void setInitialProfilePoints(int index, StsPoint pickedPoint)
    {
        if (initialProfilePoints == null)
            initialProfilePoints = StsPoint.copy(getProfilePoints());
//        System.out.println("StsVelocityProfile.setPickPoint: " + index);
        this.indexPicked = index;
        initialPickedPoint = pickedPoint;
        getSetVvsInitialProfilePoints();
    }

    public void clearInitialProfilePoints()
    {
        initialProfilePoints = null;
        initialPickedPoint = null;
        indexPicked = -99;
    }

    public void clearVvsInitialProfilePoints()
    {
        vvsInitialProfilePoints = null;
    }

    /**
     * If vvsInitialProfilePoints already set, then we will use that.
     * Otherwise, if profile points available, set them as the vvs initial profile.
     */
    public StsPoint[] getSetVvsInitialProfilePoints()
    {
        if (vvsInitialProfilePoints != null) return vvsInitialProfilePoints;
        if (profilePoints.length == 0) return null;
        vvsInitialProfilePoints = StsPoint.copy(getProfilePoints());
        return vvsInitialProfilePoints;
    }

    public void pick(GL gl, StsGLPanel glPanel)
    {
        gl.glLineWidth(2.0f);
        StsGLDraw.drawLineStrip2d(gl, StsColor.BLACK, getProfilePoints(), 2);
        gl.glLineWidth(1.0f);
    }

    public void setSelected()
    {
        System.out.println("Set this profile to the selected profile");
    }

    public StsPoint[] getProfilePoints()
    {
        if (profilePoints == null) return new StsPoint[0];
        if (lineSet != null && lineSet.getDatumShift() != 0.0)
            return correctedPoints(profilePoints);
        else
            return profilePoints;
    }

    public StsPoint[] getInitialProfilePoints()
    {
        if (initialProfilePoints == null)
            initialProfilePoints = StsPoint.copy(getProfilePoints());
        return initialProfilePoints;
    }

    public StsPoint correctedPoint(StsPoint point)
    {
        StsPoint correctedPoint = new StsPoint(point);
        correctedPoint.shiftY((float) lineSet.getDatumShift());
        return correctedPoint;
    }

    public StsPoint[] correctedPoints(StsPoint[] points)
    {
        StsPoint[] correctedPoints = StsPoint.copy(points);
        StsPoint.shiftY(correctedPoints, (float) lineSet.getDatumShift());
        return correctedPoints;
    }

    public void moveVertex(int index, StsPoint pick)
    {
        if (lineSet != null && lineSet.getDatumShift() != 0.0)
        {
            StsPoint pickCopy = new StsPoint(pick);
            pickCopy.shiftY(-(float) lineSet.getDatumShift());
            if (index < 0)
            {
                if (index == StsGather.TOP_MUTE_INDEX)
                    topMute = pickCopy;
                else if (index == StsGather.BOT_MUTE_INDEX)
                    bottomMute = pickCopy;
            }
            else
                getProfilePoints()[index] = pickCopy;
        }
    }

    public int getNProfilePoints()
    {
        return getProfilePoints().length;
    }

    public StsPoint getTopMute()
    {
        if (lineSet != null && lineSet.getDatumShift() != 0.0)
            return correctedPoint(topMute);
        else
            return topMute;
    }

    public StsPoint getBottomMute()
    {
        if (lineSet != null && lineSet.getDatumShift() != 0.0)
            return correctedPoint(bottomMute);
        else
            return bottomMute;
    }

    //	public float getStretchMute() { return stretchMute; }
    public int addProfilePoint(StsPoint pick)
    {
        if (pick == null) return -1;
        lineSet = StsPreStackLineSetClass.currentProjectPreStackLineSet;
        if (lineSet != null && lineSet.getDatumShift() != 0.0)
            pick.shiftY(-(float) lineSet.getDatumShift());
        changeType = CHANGE_POINT;
        typePicked = StsVelocityAnalysisEdit.TYPE_GATHER_EDGE;
        try
        {
            pick.v[1] = StsMath.minMax(pick.v[1], lineSet.getZMin(), lineSet.getZMax());  //make sure pick is within bounds of trace data
            if (getProfilePoints() == null || getProfilePoints().length == 0)
            {
                return addProfilePoint(pick, 0);
            }
            else
            {
                float z1 = -largeFloat;
                float z = pick.v[1];
                int nPoints = getProfilePoints().length;
                for (int n = 0; n < nPoints; n++)
                {
                    float z0 = z1;
                    z1 = getProfilePoints()[n].v[1];
                    float diff = Math.abs(z1 - z);
                    float thresh = lineSet.semblanceRangeProperties.velPickThreshold;
                    if (diff < thresh)
                    {
                        //moveVertex(n, pick); //this isn't working for right now
                        setProfilePoint(n, pick); //user clicked close to current point, so change that instead
                        typePicked = StsVelocityAnalysisEdit.TYPE_GATHER_VERTEX;
                        return n;
                    }
                    if (z >= z0 && z < z1)
                    {
                        return addProfilePoint(pick, n);
                    }
                }
                return addProfilePoint(pick, nPoints);
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsVelocityProfile.createInsertVertex() failed.", e, StsException.WARNING);
            return -1;
        }
    }

    private int addProfilePoint(StsPoint point, int arrayIndex)
    {
        setInitialProfilePoints(arrayIndex, point);
        if (debug)
            System.out.println("Adding point to profile vector of length " + getProfilePoints().length + " at index " + arrayIndex);
        setProfilePoints((StsPoint[]) StsMath.arrayInsertElementBefore(profilePoints, point, arrayIndex));
        //        fieldChanged("profilePoints", profilePoints);
        return arrayIndex;
    }

    /**
     * If this is the first pick, we already have points from an interpolated profile which have been isVisible,
     * so clear them and commit this profile.  The commit of the profile is an instanceAdd.  In addition, we are
     * going to commit the first pick with an arrayInsertCmd which adds the pick to the array of profilePoints.
     * Subsequent picks will add additional points (perhaps inserting them between existing points) to the profilePoints.
     *
     * @param typePicked  if a point is picked, this is an edgePick (point inserted between existing points) or a pointPick (point moved)
     * @param indexPicked index in array where point is to be added or inserted
     */
    public void addPickTransactionCmd(int typePicked, int indexPicked)
    {
        switch (indexPicked)
        {
            case StsGather.TOP_MUTE_INDEX:
            {
                StsChangeCmd changeCmd = new StsChangeCmd(this, topMute, "topMute", false);
                currentModel.addTransactionCmd("profile[" + getIndex() + "].topMuteChanged", changeCmd);
                break;
            }

            case StsGather.BOT_MUTE_INDEX:
            {
                StsChangeCmd changeCmd = new StsChangeCmd(this, bottomMute, "bottomMute", false);
                currentModel.addTransactionCmd("profile[" + getIndex() + "].bottomMuteChanged", changeCmd);
                break;
            }

            case StsGather.STRETCH_MUTE_INDEX:
            {
                break;
            }

            default:
                StsPoint pick = getProfilePoints()[indexPicked];
                StsDBCommand editCmd = null;
                if (typePicked == StsVelocityAnalysisEdit.TYPE_GATHER_VERTEX)
                    editCmd = new StsArrayChangeCmd(this, pick, "profilePoints", indexPicked, false);
                else
                    editCmd = new StsArrayInsertCmd(this, pick, "profilePoints", indexPicked, false);
                currentModel.addTransactionCmd("profile[" + getIndex() + "].profilePoints[" + indexPicked + "]", editCmd);
                break;
        }
        currentModel.commit();
    }

    public boolean addProfilePoint(float velocity, float y)
    {
        try
        {
            StsPoint pick = new StsPoint(velocity, y);
            setProfilePoints((StsPoint[]) StsMath.arrayAddSortedElement(profilePoints, pick, pointComparator));
            changeType = CHANGE_POINT;
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsVelocityProfile.createInsertVertex() failed.", e, StsException.WARNING);
            return false;
        }
    }

    public boolean deleteProfilePoint(StsPoint point)
    {
        int arrayIndex = StsMath.arrayGetIndex(getProfilePoints(), point);
        return deleteProfilePoint(arrayIndex);
    }

    public boolean deleteProfilePoint(int index)
    {
        setProfilePoints((StsPoint[]) StsMath.arrayDeleteElement(profilePoints, index));
        changeType = StsVelocityProfile.CHANGE_POINT;
        if (getProfilePoints().length == 0)
            delete(); // eventually does model.delete(stsObject) which creates DBCmd to delete from model
        else
        {
            StsArrayDeleteCmd cmd = new StsArrayDeleteCmd(this, "profilePoints", index, false);
            currentModel.addTransactionCmd(toString() + "." + "profilePoints", cmd);
        }
        // Since it's possible the user will individually delete all points, we need to commit each operation separately.
        // IFf we combine the deletion of all points and the profile in one transaction, then the velocityProfile will
        // have an index of -1 and writing the arrayDeleteCmd for each point will be undefined,
        // since the parent object itself is undefined (non-persistent).
        // TODO remove dysfunctional commands from the cmdList before processing
        currentModel.commit();
        return true;
    }

    public boolean hasProfilePoints()
    {
        return getProfilePoints() != null && getProfilePoints().length > 0;
    }

    public boolean clearProfilePoints()
    {
        setProfilePoints(new StsPoint[0]);
        return true;
    }

    public boolean delete()
    {
        StsPreStackLineSet lineSet = StsPreStackLineSetClass.currentProjectPreStackLineSet;
        if (lineSet == null) return false;
        if (lineSet instanceof StsPreStackLineSet3d)
            currentModel.clearDisplayTextured3dCursors(lineSet);
        else
        {
            StsPreStackLine2d line2d = (StsPreStackLine2d) ((StsPreStackLineSet2d) lineSet).getDataLine(row, col);
            if (line2d != null) line2d.textureChanged();
        }
        StsPreStackVelocityModel velocityModel = lineSet.velocityModel;
        if (velocityModel != null)
        {
            velocityModel.deleteVelocityProfile(this);
        }
        return super.delete();
    }

    public void setProfilePoints(StsPoint[] profilePoints)
    {
        this.profilePoints = profilePoints;
    }

    static class PointComparator implements Comparator
    {
        PointComparator()
        {
        }

        public int compare(Object o1, Object o2)
        {
            float z1 = ((StsPoint) o1).v[1];
            float z2 = ((StsPoint) o2).v[1];
            if (z2 > z1) return 1;
            else if (z2 < z1) return -1;
            else return 0;
        }

        public boolean equals(Object o1)
        {
            return false;
        }
    }
/*
	public void displayOnGather(StsGLPanel3d glPanel3d, float[] offsetAxisRange, StsSuperGather superGather)
                {
                    if(superGather == null) return;
                    superGather.displayNMOCurve(glPanel3d, offsetAxisRange);
                }
            */
/*
	private void drawLine(GL gl)
                {
                    gl.glBegin(GL.GL_LINE_STRIP);
                    int nPoints = profilePoints.length;
                    for (int n = 0; n < nPoints; n++)
                        gl.glVertex2fv(profilePoints[n].v, 0);
                    gl.glEnd();
                }
            */

    public void pickOnGather(StsGLPanel3d glPanel3d, StsViewGather3d viewGather)
    {
        pickOnGather(glPanel3d, (StsViewGather) viewGather);
    }

    public void pickOnGather(StsGLPanel3d glPanel3d, StsViewGather2d viewGather)
    {
        pickOnGather(glPanel3d, (StsViewGather) viewGather);
    }

    public void pickOnGather(StsGLPanel3d glPanel3d, StsViewGather viewGather)
    {
        StsSuperGather gather = viewGather.superGather;
        if (gather == null) return;
        float[] offsetAxisRange = viewGather.axisRanges[0];
        gather.pickNMOCurve(glPanel3d, offsetAxisRange);
    }

    /**
     * Pick on a gather is an offset-time point; convert to a velocity-time point; time is offset time for unflattened and tzero for flattened.
     * <p/>
     * Picking on the gather is this:
     * <p/>
     * picked offset time:                           tp = pick.getY();
     * initial velocity (before editing):            vi = initialPickedPoint.getX();
     * initial t0 (constant while bending residual): t0 = initialPickedPoint.getY();
     * offset distance:                              x = gather.getOffsetForTypeValue(pick.getX());
     * offset time at picked point:                  t = Math.sqrt(t0*t0 + x*x/(vi*vi)) + tp - t0;         equation 1
     * computed velocity for picked point:           v = (float)(Math.abs(x)/Math.sqrt(t*t - t0*t0));     equation 2
     * store new velocity in vertex:                 vertex.setX(v);
     * <p/>
     * Drawing the NMOed curve (see StsGather.displayPickNMOPoint) is this:
     * <p/>
     * velocity at picked point:                    v = vertex.point.getX();
     * offset distance at draw point:               x = drawOffsets[i];
     * initial velocity before picking:             vi = initialPickedPoint.getX();
     * offset time:                                 t = Math.sqrt(t0*t0 + x*x/(v*v));  equation 3
     * flattened offset time:                       tf = Math.sqrt(t0*t0 + x*x/(v*v)) - Math.sqrt(t0*t0 + x*x/(vi*vi) + t0;
     * <p/>
     * equations 1 and 2 should be equivalent to 3, but aren't exactly so.
     */

    public void adjustVelocityOnGather(StsPoint vertex, int index, StsPoint pick, StsSuperGather gather)
    {
        double tp = pick.getY(); // time pick on flattened or unflattened gather
        float x = Math.abs(gather.getOffsetForTypeValue(pick.getX())); // offset at pick
        double t0 = getInitialPoint(index).getY(); // tzero for NMO curve
        boolean flattened = gather.lineSet.lineSetClass.getFlatten();
        float v = 0;
        if (flattened) // compute tOffset from tzero
            tp = gather.centerGather.velocities.computeTOffset(x, tp);
        if (tp < t0 + 0.1f) return;
        v = (float) (Math.abs(x) / Math.sqrt(tp * tp - t0 * t0));
        if (debug)
        {
            System.out.println("velocity adjusted from " + vertex.getX() + " to " + v + " for t0i " + t0 + " offset " + x);
        }
        vertex.setX(v);
    }

    public StsPoint getInitialPoint(int index)
    {
        if (index >= 0)
            return initialProfilePoints[index];
        else if (index == StsGather.TOP_MUTE_INDEX)
            return topMute;
        else if (index == StsGather.BOT_MUTE_INDEX)
            return bottomMute;
        else
            return null;
    }

    public float getInitialVelocity(int index)
    {
        if (index >= 0)
            return initialProfilePoints[index].v[0];
        else if (index == StsGather.TOP_MUTE_INDEX)
            return topMute.v[0];
        else if (index == StsGather.BOT_MUTE_INDEX)
            return bottomMute.v[0];
        else
            return largeFloat;
    }

    public void adjustTzeroOnGather(StsPoint vertex, StsPoint pick)
    {
        float t0 = pick.getY();
        if (debug) System.out.println("t0 adjusted from " + vertex.getY() + " to " + t0);
        vertex.setY(t0);
        //initializeMutes();
        //		changed = true;
    }

    /**
     * linearly interpret between points in profile with constant velocity above first point and
     * extrapolated velocity below.
     *
     * @param t time
     * @return velocity at time t
     */
    public float getVelocityFromTzero(float t)
    {
        int nPoints = getProfilePoints().length;
        if (nPoints == 0) return largeFloat;
        StsPoint[] profilePoints = getProfilePoints();
        float t0 = profilePoints[0].v[1];
        float v0 = profilePoints[0].v[0];
        if (t <= t0 || nPoints == 1) return v0;
        float t1 = profilePoints[1].v[1];
        float v1 = profilePoints[1].v[0];
        int n = 2;
        while ((t > t1 || t1 == t0) && n < nPoints)
        {
            t0 = t1;
            t1 = profilePoints[n].v[1];
            v0 = v1;
            v1 = profilePoints[n].v[0];
            n++;
        }
        if (t0 == t1) return v1;
        float f = (t - t0) / (t1 - t0);
        return v0 + f * (v1 - v0);
    }

    static int nMaxIterations = 10;

    public double getVelocityFromT(float tx, float offset)
    {
        double t0 = getTZeroFromT(tx, offset, getProfilePoints());
        return offset / Math.sqrt(tx * tx - t0 * t0);
    }

    public double getInitialVelocityFromT(float tx, float offset)
    {
        double t0 = getTZeroFromT(tx, offset, getInitialProfilePoints());
        return offset / Math.sqrt(tx * tx - t0 * t0);
    }

    public double getInitialVelocityFromTzero(float z)
    {
        StsPoint[] initialProfilePoints = getInitialProfilePoints();
        if (initialProfilePoints.length == 0) return largeFloat;
        return StsMath.interpolateValue(initialProfilePoints, z, 1, 0);
    }

    private double getTZeroFromT(float tx, float offset, StsPoint[] profilePoints)
    {
        if (profilePoints.length == 0) return largeFloat;
        if (profilePoints.length == 1) return profilePoints[0].v[1];
        double v1 = profilePoints[0].v[0];
        double tz1 = profilePoints[0].v[1];
        double tx1 = Math.sqrt(tz1 * tz1 + offset * offset / (v1 * v1));
        int nPoints = profilePoints.length;
        double xsq = offset * offset;
        double txsq = tx * tx;
        double dtzdv;
        double tzFactor;
        double tzCor;
        double xFactor;
        double cor;
        for (int n = 1; n < nPoints; n++)
        {
            double v0 = v1;
            double tz0 = tz1;
            double tx0 = tx1;
            v1 = profilePoints[n].v[0];
            tz1 = profilePoints[n].v[1];
            tx1 = Math.sqrt(tz1 * tz1 + xsq / (v1 * v1));
            if (tx < tx1 || n == nPoints - 1)
            {
                dtzdv = (tz1 - tz0) / (v1 - v0);
                tzFactor = dtzdv * v0 - tz0;
                xFactor = dtzdv * dtzdv * xsq;
                double tzLast = tz0;
                double tz = tz0;
                int nIterations = 0;
                double v;
                double error;
                while (true)
                {
                    tzCor = tz + tzFactor;
                    double tzCorSq = tzCor * tzCor;
                    cor = (tz * tzCorSq * tzCor - xFactor) / ((4 * tz + tzFactor) * tzCorSq);
                    tz += cor;
                    error = Math.abs(tzLast - tz);
                    if (error < 0.01) break;
                    nIterations++;
                    if (nIterations == nMaxIterations)
                    {
                        System.out.println("getVelocityFromT failed to converge in " + nMaxIterations + " iterations. Error: " + error);
                        break;
                    }
                }
                return tz;
            }
        }
        return largeFloat;
    }

    /* old version
        private double getVelocityFromT(float tx, float offset, StsPoint[] profilePoints)
        {
            if(profilePoints.length == 0) return StsParameters.largeFloat;
            if(profilePoints.length == 1) return profilePoints[0].v[1];
            double v1 = profilePoints[0].v[0];
            double tz1 = profilePoints[0].v[1];
            double tx1 = Math.sqrt(tz1*tz1 + offset*offset/(v1*v1));
            int nPoints = profilePoints.length;
            double xsq = offset*offset;
            double txsq = tx*tx;
            for(int n = 1; n < nPoints; n++)
            {
                double v0 = v1;
                double tz0 = tz1;
                double tx0 = tx1;
                v1 = profilePoints[n].v[0];
                tz1 = profilePoints[n].v[1];
                tx1 = Math.sqrt(tz1*tz1 + xsq/(v1*v1));
                if(tx < tx1 || n == nPoints - 1)
                {
                    double f = (tx - tx0)/(tx1 - tx0);
                    double vLast = StsParameters.largeDouble;
                    int nIterations = 0;
                    double v;
                    double error;
                    while(true)
                    {
                        v = v0 + f*(v1 - v0);
                        error = Math.abs(vLast - v);
                        if(error < 0.01) break;
                        vLast = v;
                        double tz = Math.sqrt(txsq - xsq/(v*v));
                        f = (tz - tz0)/(tz1 - tz0);
                        nIterations++;
                        if(nIterations == nMaxIterations)
                        {
                            System.out.println("getVelocityFromT failed to converge in " + nMaxIterations + " iterations. Error: " + error);
                            break;
                        }
                    }
                    return v;
                }
            }
            return StsParameters.largeFloat;
        }
    */
    public double getTZeroFromT(float t, float offset)
    {
        return getTZeroFromT(t, offset, getProfilePoints());
    }

    public void checkInitializeMutes(StsPreStackLineSet lineSet)
    {
        if (topMute == null) initializeMutes(lineSet);
    }

    public void initializeMutes(StsPreStackLineSet lineSet)
    {
        StsSemblanceRangeProperties rangeProperties = lineSet.semblanceRangeProperties;
        float zMin = rangeProperties.getZMin();
        float zMax = rangeProperties.getZMax();
        if (getProfilePoints().length == 0)
        {
            topMute = new StsPoint(largeFloat, zMin);
            bottomMute = new StsPoint(largeFloat, zMax);
        }
        else
        {
            topMute = new StsPoint(getVelocityFromTzero(zMin), zMin);
            bottomMute = new StsPoint(getVelocityFromTzero(zMax), zMax);
        }
    }

    /*
public StsPoint[] getVelocities(int nCroppedSlices, float zMin, float zInc)
        {
            return StsMath.interpolatedPoints(profilePoints, nCroppedSlices, zMin, zInc, 1, 0);
        }
    */
    public int compareTo(Object other)
    {
        StsVelocityProfile otherProfile = ((StsVelocityProfile) other);
        if (otherProfile.row != row)
            return row - otherProfile.row;
        else
            return col - otherProfile.col;
    }

    public boolean sameRowCol(StsVelocityProfile otherProfile)
    {
        if (otherProfile == null) return false;
        return row == otherProfile.row && col == otherProfile.col;
    }

    public boolean addProfileOk()
    {
        return getProfilePoints().length > 2;
    }

    public String toString()
    {
        float bottomMuteTime;
        int pointsLength;


        if (bottomMute != null)
            bottomMuteTime = bottomMute.v[1];
        else
            bottomMuteTime = StsParameters.nullValue;

        if (getProfilePoints() != null)
            pointsLength = getProfilePoints().length;
        else
            pointsLength = 0;

        return "index " + getIndex() + " row " + row + " col " + col + " nValues " + pointsLength + " botMute.time " + bottomMuteTime;
    }

    static public void debugPrintProfile(Object instance, String method, StsVelocityProfile velocityProfile)
    {
        String instanceName = StsToolkit.getSimpleClassname(instance);
        System.out.print(instanceName + "." + method);
        if (velocityProfile == null)
        {
            System.out.println(" profile is null");
            return;
        }
        StsPoint[] profilePoints = velocityProfile.getProfilePoints();
        int nPoints = profilePoints.length;
        System.out.print(" Profile[" + velocityProfile.row + "][" + velocityProfile.col + "]:");
        for (int n = 0; n < nPoints; n++)
            System.out.print(" point[" + n + "]: " + profilePoints[n].v[1] + ", " + profilePoints[n].v[0]);
        System.out.println();
    }

    public static void displayTraceVelocity(StsGLPanel3d glPanel3d, StsSeismicVolume seismicVolume, float[] velocities, StsColor stsColor)
    {
        float width = 2.0f;
        GL gl = glPanel3d.getGL();
        gl.glLineWidth(width);
        float zMin = seismicVolume.getZMin();
        float zInc = seismicVolume.getZInc();
        gl.glLineStipple(1, StsGraphicParameters.dottedLine);
        gl.glEnable(GL.GL_LINE_STIPPLE);
        try
        {
            gl.glDisable(GL.GL_LIGHTING);

            stsColor.setGLColor(gl);
            float z = zMin;
            gl.glBegin(GL.GL_LINE_STRIP);
            for (int n = 0; n < velocities.length; n++, z += zInc)
            {
                gl.glVertex2f(velocities[n] / 1000, z);
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawLineStrip() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glDisable(GL.GL_LINE_STIPPLE);
            gl.glEnable(GL.GL_LIGHTING);
        }
    }

    public static void displayTraceVelocity(StsGLPanel3d glPanel3d, float[][] velocities, StsColor stsColor)
    {
        float width = 2.0f;
        GL gl = glPanel3d.getGL();
        gl.glLineWidth(width);
        gl.glLineStipple(1, StsGraphicParameters.dottedLine);
        gl.glEnable(GL.GL_LINE_STIPPLE);
        try
        {
            gl.glDisable(GL.GL_LIGHTING);

            stsColor.setGLColor(gl);
            gl.glBegin(GL.GL_LINE_STRIP);
            for (int n = 0; n < velocities.length; n++)
                gl.glVertex2fv(velocities[n], 0);

        }
        catch (Exception e)
        {
            StsException.outputException("StsGLDraw.drawLineStrip() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glDisable(GL.GL_LINE_STIPPLE);
            gl.glEnable(GL.GL_LIGHTING);
        }
    }

    /*
     public static void drawInterpolatedProfile(StsGLPanel3d glPanel3d, StsPoint[] points, StsColor stsColor)
        {
            float width = 2.0f;
            GL gl = glPanel3d.getGL();
            gl.glLineWidth(width);
            gl.glLineStipple(1, StsGraphicParameters.dottedLine);
            gl.glEnable(GL.GL_LINE_STIPPLE);
            try
            {
                gl.glDisable(GL.GL_LIGHTING);

                StsColor.setGLColor(gl, stsColor);
                gl.glBegin(GL.GL_LINE_STRIP);
                for(int n = 0; n < points.length; n++)
                    gl.glVertex2f(points[n].v[0], points[n].v[1]);
            }
            catch(Exception e)
            {
                StsException.outputException("StsGLDraw.drawLineStrip() failed.", e, StsException.WARNING);
            }
            finally
            {
                gl.glEnd();
                gl.glLineWidth(1.0f);
                gl.glDisable(GL.GL_LINE_STIPPLE);
                gl.glEnable(GL.GL_LIGHTING);
            }
        }
    */
    public static void displayNeighborProfiles(int currentRow, int currentCol, StsGLPanel3d glPanel3d,
                                               int ilineDisplayInc, int xlineDisplayInc, boolean showLabels,
                                               StsPreStackLineSet lineSet, int pickedIdx, boolean drawBackground)
    {
        StsVelocityProfile velocityProfile;

        if (ilineDisplayInc <= 0 && xlineDisplayInc <= 0) return;

        StsVelocityProfile[] neighborProfiles = lineSet.velocityModel.getNeighborProfiles(currentRow, currentCol, ilineDisplayInc, xlineDisplayInc);
        int nNeighborProfiles = neighborProfiles.length;
        for (int n = 0; n < nNeighborProfiles; n++)
            displayOnSemblance(neighborProfiles[n], glPanel3d, StsColor.GREY, false, false, showLabels, pickedIdx, drawBackground);
    }

    public static StsVelocityProfile pickNeighborProfiles(int currentRow, int currentCol, StsGLPanel3d glPanel3d, int ilineDisplayInc, int xlineDisplayInc, StsPreStackLineSet lineSet)
    {
        if (ilineDisplayInc <= 0 && xlineDisplayInc <= 0) return null;
        StsVelocityProfile[] neighborProfiles = lineSet.velocityModel.getNeighborProfiles(currentRow, currentCol, ilineDisplayInc, xlineDisplayInc);
        return (StsVelocityProfile) StsJOGLPick.pickClass3d(glPanel3d, neighborProfiles, StsJOGLPick.PICKSIZE_MEDIUM, StsJOGLPick.PICK_FIRST);
    }

    /*
        public static void displayProfile(StsGLPanel3d glPanel3d, StsSuperGather superGather, StsVelocityProfile currentVelocityProfile)
        {
            GL gl = glPanel3d.getGL();
            if(currentVelocityProfile == null) return;
            if(!superGather.centerGather.velocities.checkComputeVelocities(superGather.velocityProfile)) return;
            try
            {
                gl.glDisable(GL.GL_LIGHTING);
                gl.glLineWidth(4.0f);
                StsColor.RED.setGLColor(gl);
                gl.glBegin(GL.GL_LINE_STRIP);

                double t = superGather.seismicVolumeZMin;
                double tInc = superGather.centerGather.velocities.tInc;
                double[] vrms = superGather.centerGather.velocities.vrms;
     //           int nMin = superGather.semblanceZMinIndex;
     //           int nMax = superGather.semblanceZMaxIndex;
                for(int n = 0; n < vrms.length; n++, t += tInc)
                    gl.glVertex2d(vrms[n], t);
            }
            catch (Exception e)
            {
                StsException.outputException("StsVelocityProfile.displayIntervalVelocities() failed.", e, StsException.WARNING);
            }
            finally
            {
                gl.glEnd();
                gl.glLineWidth(1.0f);
                gl.glEnable(GL.GL_LIGHTING);
            }
        }
    */
    public static void displayIntervalVelocities(StsGLPanel3d glPanel3d, StsSuperGather superGather, StsVelocityProfile currentVelocityProfile)
    {
        GL gl = glPanel3d.getGL();
        if (currentVelocityProfile == null) return;
        if (!superGather.checkComputeVelocities()) return;
        try
        {
            gl.glDisable(GL.GL_LIGHTING);
            gl.glLineWidth(2.0f);
            StsColor.BLACK.setGLColor(gl);
            gl.glBegin(GL.GL_LINES);

            double t1 = superGather.seismicVolumeZMin;
            double tInc = superGather.centerGather.velocities.tInc;
            double[] vi = superGather.centerGather.velocities.vInt;
            //           int nMin = superGather.semblanceZMinIndex;
            //           int nMax = superGather.semblanceZMaxIndex;
            for (int n = 0; n < vi.length; n++)
            {
                double t0 = t1;
                t1 += tInc;
                gl.glVertex2d(vi[n], t0);
                gl.glVertex2d(vi[n], t1);
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsVelocityProfile.displayIntervalVelocities() failed.", e, StsException.WARNING);
        }
        finally
        {
            gl.glEnd();
            gl.glLineWidth(1.0f);
            gl.glEnable(GL.GL_LIGHTING);
        }
    }

    public static void displayOnSemblance(StsGLPanel3d glPanel3d, StsColor stsColor, boolean drawVertices, boolean active, boolean drawLabels, StsPoint[] displayPoints, int pickedIdx)
    {
        displayOnSemblance(glPanel3d, stsColor, displayPoints, drawVertices, active, drawLabels, pickedIdx, false);
    }

    public static void displayOnSemblance(StsVelocityProfile velocityProfile, StsGLPanel3d glPanel3d, StsColor stsColor, boolean drawVertices, boolean drawLabels, int pickedIdx, boolean drawBackground)
    {
        displayOnSemblance(velocityProfile, glPanel3d, stsColor, drawVertices, true, drawLabels, pickedIdx, drawBackground);
    }

    public static void displayOnSemblance(StsVelocityProfile velocityProfile, StsGLPanel3d glPanel3d, StsColor stsColor, boolean drawVertices, boolean active, boolean drawLabels, int pickedIdx, boolean drawBackground)
    {
        StsPoint[] profilePoints = velocityProfile.getProfilePoints();
        displayOnSemblance(glPanel3d, stsColor, profilePoints, drawVertices, active, drawLabels, pickedIdx, drawBackground);
    }

    public static void displayOnSemblance(StsGLPanel3d glPanel3d, StsColor stsColor, StsPoint[] profilePoints, boolean drawVertices, boolean active, boolean drawLabels, int pickedIdx, boolean drawBackground)
    {
        float width;
        if (active == true)
            width = 4.0f;
        else
            width = 2.0f;

        GL gl = glPanel3d.getGL();
        gl.glLineWidth(width);
        StsGLDraw.drawLineStrip2d(gl, StsColor.BLACK, profilePoints, 2);
        stsColor.setGLColor(gl);
        gl.glLineWidth(width / 2.0f);
        StsGLDraw.drawLineStrip2d(gl, stsColor, profilePoints, 2);
        gl.glLineWidth(1.0f);
        if (!drawVertices) return;

        drawVertices(profilePoints, gl, pickedIdx);
        if (drawLabels) drawLabels(profilePoints, gl, drawBackground);
    }

    public static void drawVertices(StsPoint[] profilePoints, GL gl, int pickedIndex)
    {
        int nPoints = profilePoints.length;
        for (int n = 0; n < nPoints; n++)
        {
            if (n == pickedIndex)
            {
                VERTEX_BORDER_COLOR.setGLColor(gl);
                StsGLDraw.drawPoint2d(profilePoints[n], gl, 12);
                VERTEX_COLOR.setGLColor(gl);
                StsGLDraw.drawPoint2d(profilePoints[n], gl, 8);
            }
            else
            {
                VERTEX_BORDER_COLOR.setGLColor(gl);
                StsGLDraw.drawPoint2d(profilePoints[n], gl, 10);
                VERTEX_COLOR.setGLColor(gl);
                StsGLDraw.drawPoint2d(profilePoints[n], gl, 6);
            }
        }
    }

    public static void drawLabels(StsPoint[] profilePoints, GL gl, boolean drawBackground)
    {
        DecimalFormat labelFormat = new DecimalFormat("####"); //decimal point not necessary SWC 6/17/09
        int nPoints = profilePoints.length;
        float shift = 0.3f;
        LABEL_COLOR.setGLColor(gl);
        for (int n = 0; n < nPoints; n++)
        {
            String text = new String(labelFormat.format(profilePoints[n].v[1]) + ", " + labelFormat.format(profilePoints[n].v[0] * 1000)); //processors prefer time/vel format SWC 6/17/09
            if (drawBackground)
                StsGLDraw.fontHelvetica12WithBackground(gl, profilePoints[n].getXYZorT(), text);
            else
                //StsGLDraw.fontHelvetica12(gl, profilePoints[n].getXYZorT(), text);
                StsGLDraw.fontOutput(gl, profilePoints[n].getX() + shift, profilePoints[n].getY(), text, com.magician.fonts.GLHelvetica18BitmapFont.getInstance(gl)); //larger text - easier to read
        }
    }

    /** sorts profile time/vel pairs by time */
    public void sort()
    {
        StsPoint.setCompareIndex(1);
        Arrays.sort(profilePoints);
    }

    /** method to get rid of points with the same time / keeps second instance */
    public void unique()
    {
        for (int i = 0; i < profilePoints.length - 1; i++)
        {
            if (Math.abs(profilePoints[i].v[1] - profilePoints[i + 1].v[1]) < lineSet.getZInc())
            {
                profilePoints = (StsPoint[]) StsMath.arrayDeleteElement(profilePoints, i);
            }
        }

    }

    /**
     * returns profile point closest in time to "time"
     * returns null if profilePoints is null
     *
     * @param time
     * @return
     */
    public StsPoint getClosestProfilePoint(double time)
    {
        if (profilePoints == null) return null;
        if (profilePoints.length == 0) return null;
        double closestDiff = Math.abs(profilePoints[0].getY() - time);
        int closestIndex = 0;
        for (int i = 1; i < profilePoints.length; i++)
        {
            double thisDiff = Math.abs(profilePoints[i].getY() - time);
            if (thisDiff < closestDiff)
            {
                closestIndex = i;
                closestDiff = thisDiff;
            }
        }
        return profilePoints[closestIndex];
    }

    public static void displayIntervalVelocityLables(StsGLPanel3d glPanel3d, StsSuperGather superGather, StsVelocityProfile currentVelocityProfile)
    {
        GL gl = glPanel3d.getGL();
        LABEL_COLOR.setGLColor(gl);
        StsPoint[] profilePoints = currentVelocityProfile.profilePoints;
        for (int n = 0; n < profilePoints.length - 1; n++)
        {
            StsPoint p0 = profilePoints[n];
            StsPoint p1 = profilePoints[n + 1];
            float t1 = p1.v[1];
            float v1 = p1.v[0];
            float t0 = p0.v[1];
            float v0 = p0.v[0];
            float tInt = (t1 + t0) / 2;
            double vInt = Math.sqrt((v1 * v1 * t1 - v0 * v0 * t0) / (t1 - t0));
            String text = new String("" + (int) Math.round(vInt * 1000));
            StsGLDraw.fontOutput(gl, v0, tInt, text, com.magician.fonts.GL9x15BitmapFont.getInstance(gl)); //small text to contrast with vels
            //StsGLDraw.fontOutput(gl, v0+shift, tInt+100, text, com.magician.fonts.GLTimesRoman24BitmapFont.getInstance(gl)); //huge text
            gl.glDrawBuffer(GL.GL_FRONT_AND_BACK);
        }
    }

    public StsVelocityProfile getIntervalVelocityProfile()
    {
        StsVelocityProfile intVelProfile = (StsVelocityProfile) this.clone();
        if (profilePoints == null || profilePoints.length == 0) return intVelProfile;
        intVelProfile.profilePoints = new StsPoint[profilePoints.length];
        intVelProfile.profilePoints[0] = new StsPoint(profilePoints[0].v[0], profilePoints[0].v[1]);
        for (int n = 0; n < profilePoints.length - 1; n++)
        {
            StsPoint p0 = profilePoints[n];
            StsPoint p1 = profilePoints[n + 1];
            float t1 = p1.v[1];
            float v1 = p1.v[0];
            float t0 = p0.v[1];
            float v0 = p0.v[0];
            double vInt = Math.sqrt((v1 * v1 * t1 - v0 * v0 * t0) / (t1 - t0));
            if (vInt == Double.NaN)
            {
                StsMessage.printMessage("NaN encountered");
            }
            intVelProfile.profilePoints[n + 1] = new StsPoint((float) vInt, t1);
        }
        // if(profilePoints.length > 1) intVelProfile.profilePoints[profilePoints.length-1] = intVelProfile.profilePoints[profilePoints.length-2];
        return intVelProfile;
    }

    public float getIntervalVelocityFromTzero(float t)
    {
        int nPoints = getProfilePoints().length;
        if (nPoints == 0) return largeFloat;
        StsPoint[] profilePoints = getProfilePoints();
        float t0 = profilePoints[0].v[1];
        float v0 = profilePoints[0].v[0];
        if (t <= t0 || nPoints == 1) return v0;
        float t1 = profilePoints[1].v[1];
        float v1 = profilePoints[1].v[0];
        int n = 2;
        while ((t > t1 || t1 == t0) && n < nPoints)
        {
            t0 = t1;
            t1 = profilePoints[n].v[1];
            v0 = v1;
            v1 = profilePoints[n].v[0];
            n++;
        }
        if (t0 == t1) return v1;
        return (float) Math.sqrt((v1 * v1 * t1 - v0 * v0 * t0) / (t1 - t0));
    }
}