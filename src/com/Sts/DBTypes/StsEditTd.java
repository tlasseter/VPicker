package com.Sts.DBTypes;

/**
 * <p>Title: jS2S development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2004</p>
 * <p>Company: S2S Systems LLC</p>
 * @author Tom Lasseter
 * @version 1.0
 */

import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.util.*;

public class StsEditTd extends StsMainObject implements StsInstance3dDisplayable
{
    public StsWell well;
    public StsPoint[] tdPoints = null;
    long time;
    public float[] adjustedTimes = null;
    transient boolean saveToDB = false;

    static StsObjectPanel objectPanel = null;
    static int pointTimeIndex = StsPoint.tIndex;
    static StsEditTd nullTdEdit = new StsEditTd("none", false);

    public StsEditTd()
    {
    }

    public StsEditTd(StsWell well, String name)
    {
        this.well = well;
        setName(name);
        computeAdjustedTimes();
    }
    public StsEditTd(StsEditTd currentEditTd)
    {
        well = currentEditTd.well;
        setNameDate();
        tdPoints = (StsPoint[])StsMath.arraycopy(currentEditTd.tdPoints);
    }

    public StsEditTd(String name, boolean persistent)
    {
        super(persistent);
        setName(name);
    }

    public StsEditTd(StsWell well)
    {
        this.well = well;
        setNameDate();

    }
	public StsEditTd(StsWell well, boolean persistent)
	{
		super(persistent);
		this.well = well;
		setNameDate();
    }

    private void setNameDate()
    {
        time = System.currentTimeMillis();
        Date date = new Date(time);
        setName(date.toString());
    }

    static public StsEditTd getNullTdEdit()
    {
        return nullTdEdit;
    }

    public void addPoint(StsPoint point)
    {
        tdPoints = (StsPoint[]) StsMath.arrayAddSortedElement(tdPoints, point);
    }

    /** One tdPoint has been adjusted in time; fit a correction cubic between this point
     *  and previous point, and this point and next point with slopes of 0 at
     *  each end of the cubic correction at each end.
     *  The equation is t = t0 + 3*dt*(m - m0)**2/dm**2 - 2*dt*(m - m0)**3/dm**3
     */
    public boolean adjustWellPath(StsWell well, StsPoint selectedPoint, StsPoint pickedPoint,
                                  StsWellViewModel wellViewModel)
    {

        for (int n = 0; n < tdPoints.length; n++)
        {
            if (tdPoints[n] == selectedPoint)
            {
                float t = selectedPoint.getT();
                float dt = pickedPoint.getT() - t;
                adjustWellPath(well, n, dt, wellViewModel);
                selectedPoint.setT(t + dt);
                setNameDate();
                return true;
            }
        }
        return false;
    }

    public void adjustWellPath(StsWell well, int selectedPointIndex, float dtMax,
                                StsWellViewModel wellViewModel)
    {
		if (well == null) return;

		StsPoint point0, point1;
        int nTdPoints = tdPoints.length;
        StsPoint[] lineVertexPoints = well.getLineVertexPoints();
        if (selectedPointIndex > 0)
        {
            point0 = tdPoints[selectedPointIndex - 1];
            point1 = tdPoints[selectedPointIndex];
            if (adjustWellPath(well, point0, point1, lineVertexPoints, 0.0f, dtMax))
            {
                saveToDB = true;
            }
        }
        if (selectedPointIndex + 1 < nTdPoints)
        {
            point0 = tdPoints[selectedPointIndex];
            point1 = tdPoints[selectedPointIndex + 1];
            if (adjustWellPath(well, point0, point1, lineVertexPoints, dtMax, 0.0f))
            {
                saveToDB = true;
            }
        }
        else
        {
            point0 = tdPoints[selectedPointIndex];
            if (shiftWellPath(well, point0, lineVertexPoints, dtMax))
            {
                saveToDB = true;
            }
        }
        if (saveToDB)
        {
            computeAdjustedTimes(lineVertexPoints);
            well.computePoints();
            well.adjustMarkerTimes();
            if (wellViewModel != null) wellViewModel.display();
            if(well.wellViewModel != null) well.wellViewModel.display();
//            fieldChanged("adjustedTimes", adjustedTimes);
//           currentModel.addMethodCmd(well, "timeAdjusted", new Object[] { adjustedTimes });
        }
    }

    public void resetTimes()
    {
        if(adjustedTimes == null) return;
		if (well == null) return;
        StsPoint[] lineVertexPoints = well.getLineVertexPoints();
        for(int n = 0; n < lineVertexPoints.length; n++)
            lineVertexPoints[n].setT(adjustedTimes[n]);

        well.computePoints();
        well.adjustMarkerTimes();
        currentModel.win3dDisplay();
    }

    private void computeAdjustedTimes()
    {
		if (well == null) return;
        computeAdjustedTimes(well.getLineVertexPoints());
    }

    private void computeAdjustedTimes(StsPoint[] lineVertexPoints)
    {
        int nLinePoints = lineVertexPoints.length;
        if (adjustedTimes == null)
        {
            adjustedTimes = new float[lineVertexPoints.length];
        }
        for (int n = 0; n < nLinePoints; n++)
        {
            adjustedTimes[n] = lineVertexPoints[n].getT();
        }
    }

    public void setAdjustedTimes(float[] times)
    {
        this.adjustedTimes = times;
    }

    public float[] getAdjustedTimes()
    {
        return adjustedTimes;
    }

    private boolean adjustWellPath(StsWell well, StsPoint point0, StsPoint point1, StsPoint[] lineVertexPoints,
                                   float t0, float t1)
    {
        float m0 = point0.getM();
        float m1 = point1.getM();
        float dt = t1 - t0;
        float dm = m1 - m0;
        float a = 3 * dt / (dm * dm);
        float b = -2 * dt / (dm * dm * dm);
        boolean adjusted = false;

        int nLinePoints = lineVertexPoints.length;
        for (int n = 0; n < nLinePoints; n++)
        {
            StsPoint linePoint = lineVertexPoints[n];
            float m = linePoint.getM();
            if (m >= m0 && m <=m1)
            {
                float dmm = m - m0;
                float dtt = t0 + (a + b * dmm) * (dmm * dmm);
                float t = linePoint.getT();
				//System.out.println("T is "+t+" changed to "+(t+dtt));
                linePoint.setT(t + dtt);
                adjusted = true;
            }
        }
        return adjusted;
    }

    private boolean shiftWellPath(StsWell well, StsPoint point0, StsPoint[] lineVertexPoints, float dt)
    {
        float m0 = point0.getM();
        boolean shifted = false;
        int nLinePoints = lineVertexPoints.length;
        for (int n = 0; n < nLinePoints; n++)
        {
            StsPoint linePoint = lineVertexPoints[n];
            float m = linePoint.getM();
            if (m > m0)
            {
                float t = linePoint.getT();
                linePoint.setT(t + dt);
                shifted = true;
            }
        }
        return shifted;
    }

    public StsPoint setSelectedPoint(int pointIndex)
    {
        return tdPoints[pointIndex];
    }

    public boolean initialize(StsModel model)
    {
		if (well != null)
			well.adjustTimes(adjustedTimes);
        return true;
    }

    public void endTransaction()
    {
        if(!saveToDB) return;
        currentModel.instanceChange(this, "tdEdit");
        saveToDB = false;
    }

    public void pick(GL gl, StsGLPanel glPanel)
    {
        if (tdPoints == null)
        {
            return;
        }
        int nPoints = tdPoints.length;
        for (int n = 0; n < nPoints; n++)
        {
            gl.glPushName(n);
            StsGLDraw.drawPoint(tdPoints[n].v, StsPoint.tIndex, StsColor.WHITE, gl, 4);
            gl.glPopName();
        }
    }

    public void display(StsGLPanel3d glPanel3d)
    {
        if (tdPoints == null)
        {
            return;
        }
        int nPoints = tdPoints.length;
        for (int n = 0; n < nPoints; n++)
        {
            StsGLDraw.drawPoint(tdPoints[n].v, pointTimeIndex, StsColor.BLACK, glPanel3d, 8,
                                StsGraphicParameters.vertexShift);
            StsGLDraw.drawPoint(tdPoints[n].v, pointTimeIndex, StsColor.WHITE, glPanel3d, 4,
                                2 * StsGraphicParameters.vertexShift);
        }
    }
}
