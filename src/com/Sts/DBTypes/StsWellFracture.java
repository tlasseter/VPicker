package com.Sts.DBTypes;

import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.Utilities.*;
import com.Sts.Utilities.Triangulation.*;

import javax.media.opengl.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Jun 4, 2008
 * Time: 8:04:29 PM
 * To change this template use File | Settings | File Templates.
 */
public class StsWellFracture extends StsMainObject implements StsTreeObjectI
{
    private StsColor stsColor;
    private StsWell well;
    double[][] points;
    int[][] triangleIndexes;
    double[] normal;
    private double zMin;
    private double zMax;
    private double dipAngle;
    private double azimuth;

    static final double[] vertical = new double[]{0.0, 0.0, -1.0};

    static public StsFieldBean[] displayFields = null;
    static StsDoubleFieldBean zMinBean;
    static StsDoubleFieldBean zMaxBean;
    static StsDoubleFieldBean dipAngleBean;
    static StsDoubleFieldBean azimuthBean;

    static protected StsObjectPanel objectPanel = null;

    public StsWellFracture()
    {
    }

    private StsWellFracture(StsLstSqPointsPlane plane, String name, StsColor color, StsWell well)
    {
        super(false);
        setName(name);
        stsColor = color;
        this.well = well;
        triangulate(plane);
        computeRotatedPoints();
        addToModel();
    }

    static public StsWellFracture constructor(StsLstSqPointsPlane plane, String name, StsColor color, StsWell well)
    {
        try
        {
            return new StsWellFracture(plane, name, color, well);
        }
        catch (Exception e)
        {
            StsException.systemError(StsWellFracture.class, "constructor");
            return null;
        }
    }

    private void triangulate(StsLstSqPointsPlane plane)
    {
        normal = plane.getNormal();
        double[] axis = StsMath.cross(normal, vertical);
        setAzimuth(StsMath.atan2d(axis[1], axis[0]));
        // StsToolkit.print(axis, "Axis: ");
        double axisLen = StsMath.length(axis);
        float angleRad = (float) Math.asin(axisLen);
        float angle = (float) StsMath.DEGperRAD * angleRad;
        // System.out.println("Rotation angle: " + angle);
        dipAngle = Math.abs(90 - angle);
        StsQnion q = new StsQnion(axis, angleRad);
        points = plane.getPoints();
        int nPoints = points.length;
        double[][] rotatedPoints = new double[nPoints][];
        for (int n = 0; n < nPoints; n++)
            rotatedPoints[n] = q.leftRotateVec(points[n]);
        StsDelaunayTriangulation triangulation = new StsDelaunayTriangulation(rotatedPoints);
        triangleIndexes = triangulation.triangleIndexes;
        computeZRange();
    }

    private void computeRotatedPoints()
    {
        StsProject project = currentModel.getProject();
        for (int n = 0; n < points.length; n++)
        {
            float[] xy = project.getRotatedRelativeXYFromUnrotatedRelativeXY((float)points[n][0], (float)points[n][1]);
            points[n][0] = xy[0];
            points[n][1] = xy[1];
        }
    }

    private void computeZRange()
    {
        if (points == null) return;
        zMin = points[0][2];
        zMax = zMin;

        for (int n = 1; n < points.length; n++)
        {
            double z = points[n][2];
            if (z < zMin) zMin = z;
            else if (z > zMax) zMax = z;
        }
    }

    public void display(StsGLPanel glPanel)
    {
        if (points == null || triangleIndexes == null) return;
        GL gl = glPanel.getGL();
        try
        {
            stsColor.setGLColor(gl);
            int nTriangles = triangleIndexes.length;
            gl.glBegin(GL.GL_TRIANGLES);
            for (int n = 0; n < nTriangles; n++)
            {
                gl.glNormal3dv(normal, 0);
                int[] indexes = triangleIndexes[n];
                for (int i = 0; i < 3; i++)
                {
                    double[] point = points[indexes[i]];
                    gl.glVertex3dv(point, 0);
                }
            }
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "display.", e);
        }
        finally
        {
            gl.glEnd();
        }
    }

    public boolean initialize(StsModel model)
    {
        return true;
    }

    public StsFieldBean[] getDisplayFields()
    {
        try
        {
            if (displayFields == null)
            {
                displayFields = new StsFieldBean[]
                    {
                        zMinBean = new StsDoubleFieldBean(StsWellFracture.class, "zMin", true, "Depth Min"),
                        zMaxBean = new StsDoubleFieldBean(StsWellFracture.class, "zMax", true, "Depth Max"),
                        dipAngleBean = new StsDoubleFieldBean(StsWellFracture.class, "dipAngle", true, "Dip Angle"),
                        azimuthBean = new StsDoubleFieldBean(StsWellFracture.class, "azimuth", true, "Azimuth")
                    };
            }
            return displayFields;
        }
        catch (Exception e)
        {
            StsException.outputException("StsSensor.getDisplayFields() failed.", e, StsException.WARNING);
            return null;
        }
    }

    public StsFieldBean[] getPropertyFields()
    {
        return null;
    }

    public Object[] getChildren()
    {
        return new Object[0];
    }

    public boolean anyDependencies()
    {
        return false;
    }

    public StsObjectPanel getObjectPanel()
    {
        if (objectPanel == null)
        {
            objectPanel = StsObjectPanel.constructor(this, true);
        }
        return objectPanel;
    }

    public void treeObjectSelected()
    {
    }

    public double getZMin()
    {
        return zMin;
    }

    public void setZMin(double zMin)
    {
        this.zMin = zMin;
    }

    public double getZMax()
    {
        return zMax;
    }

    public void setZMax(double zMax)
    {
        this.zMax = zMax;
    }

    public double getDipAngle()
    {
        return dipAngle;
    }

    public void setDipAngle(double dipAngle)
    {
        this.dipAngle = dipAngle;
    }

    public double getAzimuth()
    {
        return azimuth;
    }

    public void setAzimuth(double azimuth)
    {
        this.azimuth = azimuth;
    }
}
