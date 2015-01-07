//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System


package com.Sts.DBTypes;

import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.util.*;

public class StsLine extends StsMainTimeObject implements StsTreeObjectI, StsInstance3dDisplayable, StsViewSelectable, StsTimeDepthDisplayable
{
    /** Section this line is on (if any) */
    protected StsSection onSection = null;
    /** Sections connected to this line */
    protected StsObjectRefList connectedSections = null;
    /** For a faultSection: RIGHT or LEFT */
    protected int sectionSide = 0;
    protected StsColor stsColor = new StsColor(StsColor.RED);
    protected boolean isVertical = false;
    /** StsSurfaceVertex's which define line geometry: unrotated coor system */
    protected StsObjectRefList lineVertices = null;
    /** StsSurfaceVertex at intersection of well, surface, block: rotated coor system */
    protected StsObjectRefList sectionEdgeVertices = null;
    /** Indicates lineVertices are in rotated coordinate system */
    private boolean isVerticesRotated = true;

    protected double xOrigin, yOrigin;
    protected boolean drawZones = false;
    protected boolean drawLabels = false;
    /** Indicates whether this well has depth, time, or time and depth coordinates */
    protected byte zDomainSupported = TD_DEPTH;
    /** original domain this object was built in. If and when velocity model is rebuilt, this domain is used as the coordinate source. */
    private byte zDomainOriginal = StsParameters.TD_DEPTH;

    transient StsObjectList surfaceEdgeVertices = null;

    transient protected StsList zones = null;
    transient protected StsPoint[] rotatedPoints = null; /* array of segmented points in rotated coor system */
    transient protected int iTop, iBot;
    transient protected boolean drawVertices = false;
    transient protected StsSurfaceVertex selectedVertex = null;

    /** In Bezier curve fitting, max distance between generated points */
    static float maxSplineDistance = 100.0f;

    static StsObjectPanel objectPanel = null;

    static protected boolean checkExcursion = false;

    static final float nullValue = StsParameters.nullValue;

    protected boolean highlighted = false;
    /** draw the well highlighted */

    transient static public StsHighlightedList highlightedList = new StsHighlightedList(4, 2);
    transient static protected StsLine currentSectionLine = null;

    transient protected StsBoundingBox boundingBox = null;

    // display fields: renamed to pseudoDisplayFields because parent class has static displayFields also
    static public StsFieldBean[] pseudoDisplayFields = null;

    public transient boolean initialized = false;
    /** this line has been initialized */
    public transient boolean sectionInitialized = false;
    /** true if this line has been projected to section */

    public static final int RIGHT = StsParameters.RIGHT;
    public static final int LEFT = StsParameters.LEFT;

    static public final int NONE = StsParameters.NONE;
    static public final int MINUS = StsParameters.MINUS;
    static public final int PLUS = StsParameters.PLUS;

    static public final byte TD_DEPTH = StsParameters.TD_DEPTH;
    static public final byte TD_TIME = StsParameters.TD_TIME;
    static public final byte TD_APPROX_DEPTH = StsParameters.TD_APPROX_DEPTH;
    static public final byte TD_TIME_DEPTH = StsParameters.TD_TIME_DEPTH;
    static public final byte TD_APPROX_DEPTH_AND_DEPTH = StsParameters.TD_APPROX_DEPTH_AND_DEPTH;
    static public final byte TD_NONE = StsParameters.TD_NONE;

    static public final long serialVersionUID = 1L;

    /** default constructor */
    public StsLine() // throws StsException
    {
        /*
            if(currentModel == null) return;
    //        setName(new String("PseudoLine-" + index));
            connectedSections = StsObjectRefList.constructor(2, 2, "connectedSections", this);
            lineVertices = StsObjectRefList.constructor(2, 2, "lineVertices", this);
            surfaceVertices = StsObjectRefList.constructor(2, 2, "surfaceVertices", this);
            StsProject project = currentModel.getProject();
            stsColor = new StsColor(Color.RED);
            xOrigin = project.getXOrigin();
            yOrigin = project.getYOrigin();
        */
    }

    public StsLine(String name, boolean persistent)
    {
        this(persistent);
        setName(name);
    }

    public StsLine(boolean persistent)
    {
        super(persistent);
    }

    static public void initColors()
    {
        StsColorComboBoxFieldBean colorComboBox;
        StsSpectrum spectrum = currentModel.getSpectrum("Basic");
        if (pseudoDisplayFields == null)
            return;
        colorComboBox = (StsColorComboBoxFieldBean) StsFieldBean.getBeanWithFieldName(pseudoDisplayFields, "stsColor");
        colorComboBox.setListItems(spectrum);
    }

    static public StsLine buildLine()
    {
        StsLine line = new StsLine(false);
        line.setZDomainOriginal(currentModel.getProject().getZDomain());
        line.initializeVertices();
        return line;
    }

    public void initializeVertices()
    {
        lineVertices = StsObjectRefList.constructor(2, 2, "lineVertices", this);
        sectionEdgeVertices = StsObjectRefList.constructor(2, 2, "sectionEdgeVertices", this);
        numberOfElements = lineVertices.getSize();
    }

    public String getName()
    {
        if (name != null)
            return name;
        else
            return "Line-" + getIndex();
    }

    public void setDrawZones(boolean state)
    {
        drawZones = state;
    }

    public boolean getDrawZones()
    {
        return drawZones;
    }

    public void setDrawLabels(boolean state)
    {
        drawLabels = state;
    }

    public boolean getDrawLabels()
    {
        return drawLabels;
    }

    public void setIsVertical(boolean isVertical)
    {
        this.isVertical = isVertical;
    }

    public boolean getIsVertical()
    {
        return isVertical;
    }

    public StsSection getOnSection()
    {
        return onSection;
    }

    public int getSectionSide()
    {
        return sectionSide;
    }

    public StsPoint[] getRotatedPoints() { return rotatedPoints; }

    /* public void setRotatedPoints(StsPoint[] rotatedPoints) { this.rotatedPoints = rotatedPoints; } */

    public StsPoint[] getExportPoints() { return getLineVertexPoints(); }

    public boolean normalShift(float offset)
    {
        // compute scale factors in X & Y
        float dist = StsMath.distance(rotatedPoints[0].getXYZ(), rotatedPoints[rotatedPoints.length - 1].getXYZ(), 2);
        float[] normal = StsMath.horizontalNormal(rotatedPoints[0].getPointXYZ(), rotatedPoints[rotatedPoints.length - 1].getXYZ(), 1);
        float xScale = normal[0] / dist;
        float yScale = normal[1] / dist;
        // compute shifts normal to the azimuth

        this.xOrigin = xOrigin + xScale * offset;
        this.yOrigin = yOrigin + yScale * offset;

        computePoints();
        return true;
    }

    public int getTopIndex()
    {
        return iTop;
    }

    public int getBotIndex()
    {
        return iBot;
    }

    public StsObjectRefList getSectionEdgeVertices()
    {
        return sectionEdgeVertices;
    }

    public StsObjectRefList getLineVertices()
    {
        return lineVertices;
    }

    public StsSurfaceVertex getTopLineVertex()
    {
        return (StsSurfaceVertex) lineVertices.getFirst();
    }

    public StsSurfaceVertex getBotLineVertex()
    {
        return (StsSurfaceVertex) lineVertices.getLast();
    }

    /** returns the top line vertex point.  This is a unrotated coordinate point relative to the
     *  origin point.  Don't use as a rotated local coordinate, because it isn't.  Along with the origin,
     *  it allows you to recover the absolute global coordinates of this point.
     * @return unrotated coordinate offset from well origin (typically Kelly bushing).
     */
    public StsPoint getTopPoint()
    {
        return getTopLineVertex().getPoint();
    }


    /** returns the bot line vertex point.  This is a unrotated coordinate point relative to the
     *  origin point.  Don't use as a rotated local coordinate, because it isn't.  Along with the origin,
     *  it allows you to recover the absolute global coordinates of this point.
     * @return unrotated coordinate offset from well origin (typically Kelly bushing).
     */
    public StsPoint getBotPoint()
    {
        return getBotLineVertex().getPoint();
    }

    public StsPoint getTopRotatedPoint()
    {
        return rotatedPoints[0];
    }

    public StsPoint getBotRotatedPoint()
    {
        return rotatedPoints[rotatedPoints.length - 1];
    }

    public StsPoint getBotVector()
    {
        if (rotatedPoints == null) return null;
        StsPoint lastPoint = rotatedPoints[rotatedPoints.length - 1];
        StsPoint prevPoint = rotatedPoints[rotatedPoints.length - 2];
        return StsPoint.subPointsStatic(lastPoint, prevPoint);
    }

    public float getTopZ()
    {
        return getTopPoint().getZorT();
    }

    public float getBotZ()
    {
        return getBotPoint().getZorT();
    }

    public byte getZDomainSupported()
    {
        return zDomainSupported;
    }

    public StsPoint[] getUnrotatedPoints() { return getLineVertexPoints(); }

    public StsPoint[] getLineVertexPoints()
    {
        if (lineVertices == null)
        {
            return null;
        }
        int nLineVertices = lineVertices.getSize();
        StsPoint[] lineVertexPoints = new StsPoint[nLineVertices];
        for (int n = 0; n < nLineVertices; n++)
        {
            lineVertexPoints[n] = ((StsSurfaceVertex) lineVertices.getElement(n)).getPoint();
        }
        return lineVertexPoints;
    }

    public void deleteVertexNearestPoint(StsPoint point)
    {
        if (lineVertices == null)
        {
            return;
        }
        int nLineVertices = lineVertices.getSize();
        StsPoint[] lineVertexPoints = new StsPoint[nLineVertices];
        for (int n = 0; n < nLineVertices; n++)
        {
            lineVertexPoints[n] = ((StsSurfaceVertex) lineVertices.getElement(n)).getPoint();
            if (point.getZorT() < lineVertexPoints[n].getZorT())
            {
                if (n == 0)
                {
                    deleteSectionVertex((StsSurfaceVertex) lineVertices.getElement(n));
                }
                else
                {
                    if ((lineVertexPoints[n].getZorT() - point.getZorT()) >
                        (point.getZorT() - lineVertexPoints[n - 1].getZorT()))
                    {
                        deleteSectionVertex((StsSurfaceVertex) lineVertices.getElement(n - 1));
                    }
                    else
                    {
                        deleteSectionVertex((StsSurfaceVertex) lineVertices.getElement(n));
                    }
                }
            }
        }
        return;
    }

    public void setHighlighted(boolean val)
    {
        highlighted = val;
        dbFieldChanged("highlighted", highlighted);
        currentModel.win3dDisplayAll();
    }

    public boolean getHighlighted() { return highlighted; }

    public void setOnSection(StsSection onSection)
    {
        this.onSection = onSection;
        // reinitialize=true: on db reload reinitialize this as section has changed
        dbFieldChanged("onSection", onSection, true);
    }

    public StsSection[] getAllSections()
    {
        int nTotalSections = 0;
        if(onSection != null)
            nTotalSections = 1;
        int nConnectedSections = 0;
        if(connectedSections != null)
        {
            nConnectedSections = connectedSections.getSize();
            nTotalSections += nConnectedSections;
        }
        StsSection[] sections = new StsSection[nTotalSections];
        int n = 0;
        if(onSection != null)
            sections[n++] = onSection;
        if(connectedSections == null) return sections;
        for(int i = 0; i < nConnectedSections; i++)
            sections[n++] = (StsSection)connectedSections.getElement(i);
        return sections;
    }

    public StsSection[] getConnectedSections()
    {
        if(connectedSections == null)
            return new StsSection[0];
        else
            return (StsSection[]) connectedSections.getCastList(StsSection.class);
    }

    public void setInitialized(boolean initialized)
    {
        this.initialized = initialized;
    }

    public boolean getInitialized()
    {
        return initialized;
    }

    public void setDrawVertices(boolean draw)
    {
        drawVertices = draw;
    }

    public boolean getDrawVertices()
    {
        return drawVertices;
    }

    public void setSelectedVertex(StsSurfaceVertex v)
    {
        selectedVertex = v;
    }

    public StsSurfaceVertex getSelectedVertex()
    {
        return selectedVertex;
    }

    public StsColor getStsColor()
    {
        return stsColor;
    }

    public void setStsColor(StsColor color)
    {
        if (stsColor == color) return;
        stsColor = color;
        dbFieldChanged("stsColor", color);
    }

    public double getXOrigin()
    {
        return xOrigin;
    }

    public void setXOrigin(double xOrigin)
    {
        this.xOrigin = xOrigin;
    }

    public double getYOrigin()
    {
        return yOrigin;
    }

    public void setYOrigin(double yOrigin)
    {
        this.yOrigin = yOrigin;
    }

    /*
       public StsSection[] getAssociatedSections()
       {
           StsSection[] associatedSections = new StsSection[0];
    if(section != null) associatedSections = (StsSection[])StsMath.arrayAddElement(associatedSections, section);
           if(connectedSections == null) return associatedSections;

           for(int n = 0; n < connectedSections.length; n++)
           {
               StsSection connectedSection = (StsSection)connectedSections.getElement(n);
               associatedSections = (StsSection[])StsMath.arrayAddElement(associatedSections, connectedSection);
           }
           return associatedSections;
       }
    */
    // there can be a maximum of two sections connections to this well/fault
    // or one connected section if well is on a section
    public boolean isFullyConnected()
    {
        int nConnectedSections = getNConnectedSections();
        if (onSection != null)
        {
            if (nConnectedSections == 1) return true;
            StsException.systemError(this, "isFullyConnected", "Line " + toSectionString() +
                " not properly connected. Is on a section and connected to " + nConnectedSections + " other sections.");
            return false;
        }
        return nConnectedSections == 2;
    }

    public int getNConnectedSections()
    {
        if (connectedSections == null) return 0;
        return connectedSections.getSize();
    }

    public boolean isFullyConnected(StsModelSurface surface)
    {
        if (connectedSections == null)
        {
            return false;
        }
        int nConnectedSections = connectedSections.getSize();
        if (nConnectedSections == 0)
        {
            return false;
        }
        for (int n = 0; n < nConnectedSections; n++)
        {
            StsSection section = (StsSection) connectedSections.getElement(n);
            if (!section.hasSurface(surface))
            {
                return false;
            }
        }
        return true;
    }

    public int getMaxConnections()
    {
        if (onSection != null)
        {
            return 1;
        }
        else
        {
            return 2;
        }
    }

    public boolean addConnectedSection(StsSection section)
    {
        int nTotalSections = 0;
        if (this.onSection != null)
            nTotalSections = 1;
        if (connectedSections == null)
            connectedSections = StsObjectRefList.constructor(2, 2, "connectedSections", this);
        else
        {
            int nConnectedSections = connectedSections.getSize();
            nTotalSections += nConnectedSections;
        }
        if (nTotalSections >= 2)
        {
            StsMessageFiles.errorMessage("Two sections already connected at this StsLine: " + toSectionString() + ". Cannot add another.");
            return false;
        }
        return connectedSections.add(section);
    }

    public void setSectionLine()
    {
        StsSection connectedSection = getOnlyConnectedSection();
        connectedSection.setLine(this);
    }

    public boolean deleteConnectedSection(StsSection section)
    {
        if (connectedSections == null)
        {
            return false;
        }
        return connectedSections.delete(section);
    }

    public StsSection getOnlyConnectedSection()
    {
        if (connectedSections == null || connectedSections.getSize() > 1)
        {
            return null;
        }
        else
        {
            return (StsSection) connectedSections.getElement(0);
        }
    }

    public String getLabel()
    {
        return getName() + " ";
    }

    public String lineOnSectionLabel()
    {
        if (onSection == null)
        {
            return getLabel();
        }
        else
        {
            return new String(getLabel() + onSection.getLabel() + " " +
                StsParameters.sideLabel(sectionSide) + " ");
        }
    }

    public int getLineSectionEnd()
    {
        int nConnectedSections = connectedSections.getSize();
        if (nConnectedSections <= 0)
        {
            return NONE;
        }
        StsSection firstConnectedSection = (StsSection) connectedSections.getElement(0);
        if (firstConnectedSection == null)
        {
            return NONE;
        }
        return firstConnectedSection.getSidePosition(this);
    }

    // If well is on an auxilarySection or connected to one: it is a dyingFault.
    public boolean isDyingFault()
    {
        if (onSection != null)
        {
            return onSection.isAuxiliary();
        }

        int nConnectedSections = connectedSections == null ? 0 : connectedSections.getSize();
        for (int n = 0; n < nConnectedSections; n++)
        {
            StsSection connectedSection = (StsSection) connectedSections.getElement(n);
            if (connectedSection.isAuxiliary())
            {
                return true;
            }
        }
        return false;
    }

    public boolean isFault()
    {
        return type == StsParameters.FAULT;
    }

    public void setSectionSide(int side)
    {
        sectionSide = side;
        dbFieldChanged("sectionSide", side);
    }

    public StsPoint getXYZPointAtZorT(float z, boolean extrapolate)
    {
        return getXYZPointAtZorT(z, extrapolate, isDepth);
    }

    /**
     * Line vertices have 5 coordinates (x, y, z, mDepth, time). Return a new Point
     * at z with 3 coordinates. z is either z or time, depending on isDepth.
     *
     * @param z           float
     * @param extrapolate boolean
     * @return StsPoint
     */
    public StsPoint getXYZPointAtZorT(float z, boolean extrapolate, boolean isDepth)
    {
        StsPoint point;
        int i = 0;
        float z0 = 0.0f, z1 = 0.0f;

        if (rotatedPoints == null)
        {
            if(!computePoints()) return null;
        }

        int nPoints = rotatedPoints.length;
        if (nPoints < 2)
        {
            point = rotatedPoints[0];
            return new StsPoint(point.getX(), point.getY(), z, point.getM());
        }

        z0 = rotatedPoints[0].getZorT(isDepth);
        if (z < z0)
        {
            if (extrapolate)
            {
                z1 = rotatedPoints[1].getZorT(isDepth);
                i = 1;
            }
            else
            {
                point = rotatedPoints[0];
                return new StsPoint(point.getX(), point.getY(), z, point.getM());
            }
        }
        else if ((z1 = rotatedPoints[nPoints - 1].getZorT(isDepth)) < z)
        {
            if (extrapolate)
            {
                z0 = rotatedPoints[nPoints - 2].getZorT(isDepth);
                i = nPoints - 1;

            }
            else
            {
                point = rotatedPoints[nPoints - 1];
                return new StsPoint(point.getX(), point.getY(), z1);
            }
        }
        else
        {
            z1 = rotatedPoints[0].getZorT(isDepth);
            for (i = 1; i < nPoints; i++)
            {
                z0 = z1;
                z1 = rotatedPoints[i].getZorT(isDepth);
                if (z >= z0 && z <= z1)
                {
                    break;
                }
            }
        }

        float ff = (z - z0) / (z1 - z0);
        //        point = StsPoint.staticInterpolatePoints(points[i - 1], points[i], ff);
        point = new StsPoint(4);
        point.interpolatePoints(rotatedPoints[i - 1], rotatedPoints[i], ff);
        point.setZ(z);
        return point;
    }

    public StsPoint getPointAtZorT(float z, boolean extrapolate)
    {
        return getPointAtZorT(z, extrapolate, isDepth);
    }

/**
 * Line vertices have 5 coordinates (x, y, z, mDepth, time). Return a new Point
      * at z with 3 coordinates. z is either z or time, depending on isDepth.
      *
      * @param z           float
      * @param extrapolate boolean
      * @return StsPoint
      */
     public StsPoint getPointAtZorT(float z, boolean extrapolate, boolean isDepth)
     {
         StsPoint point;
         int i = 0;
         float z0 = 0.0f, z1 = 0.0f;

         if (rotatedPoints == null)
         {
             if(!computePoints()) return null;
         }

         int nPoints = rotatedPoints.length;
         if (nPoints < 2)
         {
             point = rotatedPoints[0];
             return new StsPoint(point);
         }

         z0 = rotatedPoints[0].getZorT(isDepth);
         if (z < z0)
         {
             if (extrapolate)
             {
                 z1 = rotatedPoints[1].getZorT(isDepth);
                 i = 1;
             }
             else
             {
                 point = rotatedPoints[0];
                 return new StsPoint(point);
             }
         }
         else if ((z1 = rotatedPoints[nPoints - 1].getZorT(isDepth)) < z)
         {
             if (extrapolate)
             {
                 z0 = rotatedPoints[nPoints - 2].getZorT(isDepth);
                 i = nPoints - 1;

             }
             else
             {
                 point = rotatedPoints[nPoints - 1];
                 return new StsPoint(point.getX(), point.getY(), z1);
             }
         }
         else
         {
             z1 = rotatedPoints[0].getZorT(isDepth);
             for (i = 1; i < nPoints; i++)
             {
                 z0 = z1;
                 z1 = rotatedPoints[i].getZorT(isDepth);
                 if (z >= z0 && z <= z1)
                 {
                     break;
                 }
             }
         }

         float ff = (z - z0) / (z1 - z0);
         point = StsPoint.staticInterpolatePoints(rotatedPoints[i - 1], rotatedPoints[i], ff);
         //point = new StsPoint(4);
         //point.interpolatePoints(rotatedPoints[i - 1], rotatedPoints[i], ff);
         //point.setZ(z);
         return point;
     }

    public StsPoint getPointAtZ(float z, boolean extrapolate)
    {
        return getXYZPointAtZorT(z, extrapolate, true);
    }

    /*
        public StsPoint getPointAtZ(float z, boolean extrapolate)
        {
            StsPoint point;
            int i = 0;
            float z0 = 0.0f, z1 = 0.0f;

            if (points == null)
            {
                return null;
            }

            int nPoints = points.length;
            if (isVertical || nPoints < 2)
            {
                return new StsPoint(points[0].getX(), points[0].getY(), z);
            }

            if (z < points[0].getZ(isDepth))
            {
                if (extrapolate)
                {
                    z0 = points[0].getZ(isDepth);
                    z1 = points[1].getZ(isDepth);
                    i = 1;
                }
                else
                {
                    return new StsPoint(points[0].getX(), points[0].getY(), points[0].getZ(isDepth));
                }
            }
            else if (z > points[nPoints - 1].getZ(isDepth))
            {
                if (extrapolate)
                {
                    z0 = points[nPoints - 2].getZ(isDepth);
                    z1 = points[nPoints - 1].getZ(isDepth);
                    i = nPoints - 1;

                }
                else
                {
                    return new StsPoint(points[nPoints - 1].getX(), points[nPoints - 1].getY(),
                                        points[nPoints - 1].getZ(isDepth));
                }
            }
            else
            {
                z1 = points[0].getZ(isDepth);
                for (i = 1; i < nPoints; i++)
                {
                    z0 = z1;
                    z1 = points[i].getZ(isDepth);
                    if (z >= z0 && z <= z1)
                    {
                        break;
                    }
                }
            }

            float ff = (z - z0) / (z1 - z0);
            point = StsPoint.staticInterpolatePoints(points[i - 1], points[i], ff);
            return new StsPoint(point.getX(), point.getY(), z);
        }
    */
    public StsPoint getPointAtMDepth(float m, boolean extrapolate)
    {
        StsPoint[] points = getRotatedPoints();
        if(points == null) return null;
        return StsMath.interpolatePoint(m, points, 3, extrapolate);
    }

    /** if an StsLine has no measured depths, compute them by linear segments between line vertices */
    public void computeMDepths()
    {
        int nPoints = lineVertices.getSize();
        StsPoint point2 = ((StsSurfaceVertex) lineVertices.getFirst()).getPoint();
        float mdepth = 0.0f;
        point2.setM(0.0f);
        for (int n = 1; n < nPoints; n++)
        {
            StsPoint point1 = point2;
            point2 = ((StsSurfaceVertex) lineVertices.getElement(n)).getPoint();
            mdepth += StsMath.distance(point1, point2);
            point2.setM(mdepth);
        }
    }

    public void checkSortLineVertices()
    {
        StsSurfaceVertex[] surfaceVertices = (StsSurfaceVertex[])lineVertices.getCastList();
        Arrays.sort(surfaceVertices);
        lineVertices.getList().list = surfaceVertices;
    }

    public void removeClosePoints()
    {
        int nPoints = lineVertices.getSize();
        float z = ((StsSurfaceVertex)lineVertices.getLast()).getPoint().getZ();
        for(int n = nPoints-2; n >= 0; n--)
        {
            float prevZ = z;
            z = ((StsSurfaceVertex)lineVertices.getElement(n)).getPoint().getZ();
            if(Math.abs(z-prevZ) < 1.0f)
                lineVertices.delete(n);
        }
    }

    public void checkIsVertical()
    {
        int nPoints = lineVertices.getSize();
        StsPoint topPoint = ((StsSurfaceVertex)lineVertices.getFirst()).getPoint();
        for(int n = 1; n < nPoints; n++)
        {
            StsPoint point = ((StsSurfaceVertex)lineVertices.getElement(n)).getPoint();
            if(!point.sameXY(topPoint)) return;
        }
        isVertical = true;
    }

    public StsPoint getSlopeAtMDepth(float m)
    {
        return getSlopeAtMDepth(m, rotatedPoints);
    }

    //TODO this could be made more efficient by taking advantage of fact that mdepths and rotatedPoints are monotonically increasing in depth
    public StsPoint[] getSlopesAtMDepths(float[] mdepths)
    {
        int nMdepthPoints = mdepths.length;
        StsPoint[] slopes = new StsPoint[nMdepthPoints];
        for (int n = 0; n < nMdepthPoints; n++)
            slopes[n] = getSlopeAtMDepth(mdepths[n]);
        return slopes;
    }

    //TODO this could be made more efficient by taking advantage of fact that mdepths and rotatedPoints are monotonically increasing in depth
    public StsPoint[] getPointsAtMDepths(float[] mdepths)
    {
        int nMdepthPoints = mdepths.length;
        StsPoint[] slopes = new StsPoint[nMdepthPoints];
        for (int n = 0; n < nMdepthPoints; n++)
            slopes[n] = getPointAtMDepth(mdepths[n], true);
        return slopes;
    }


    public StsPoint getPointAtMDepth(float m, StsPoint[] points, boolean extrapolate)
    {
        // 3 is measured depth index in point.v array
        return StsMath.interpolatePoint(m, points, 3, extrapolate);
    }

    public StsPoint getSlopeAtMDepth(float m, StsPoint[] points)
    {
        StsPoint slope;
        int i = 0;
        float f;

        if (points == null)
        {
            return null;
        }

        int nPoints = points.length;
        if (isVertical || nPoints < 2)
        {
            slope = new StsPoint(5);
            slope.v[2] = 1.0f;
        }
        else
        {
            int index = StsMath.arrayIndexBelow(m, points, 3);
            index = StsMath.minMax(index, 0, points.length - 2);
            slope = StsPoint.subPointsStatic(points[index + 1], points[index]);
            slope.normalizeXYZ();
        }
        return slope;
    }

    /** returns a point at Z with all coordinates: x, y, z, mdepth, time */
    public StsPoint getLinePointAtZ(float z, boolean extrapolate)
    {
        StsPoint point;
        int i = 0;
        float z0 = 0.0f, z1 = 0.0f;

        if (rotatedPoints == null)
        {
            return null;
        }

        int nPoints = rotatedPoints.length;
        if (isVertical || nPoints < 2)
        {
            point = rotatedPoints[0].copy();
            //            point = new StsPoint(points[0]);
            point.setZorT(z, isDepth);
            return point;
        }

        if (z < rotatedPoints[0].getZorT())
        {
            if (extrapolate)
            {
                z0 = rotatedPoints[0].getZorT();
                z1 = rotatedPoints[1].getZorT();
                i = 1;
            }
            else
            {
                return new StsPoint(rotatedPoints[0]);
            }
        }
        else if (z > rotatedPoints[nPoints - 1].getZorT())
        {
            if (extrapolate)
            {
                z0 = rotatedPoints[nPoints - 2].getZorT();
                z1 = rotatedPoints[nPoints - 1].getZorT();
                i = nPoints - 1;

            }
            else
            {
                return new StsPoint(rotatedPoints[nPoints - 1]);
            }
        }
        else
        {
            z1 = rotatedPoints[0].getZorT();
            for (i = 1; i < nPoints; i++)
            {
                z0 = z1;
                z1 = rotatedPoints[i].getZorT();
                if (z >= z0 && z <= z1)
                {
                    break;
                }
            }
        }

        float ff = (z - z0) / (z1 - z0);
        point = StsPoint.staticInterpolatePoints(rotatedPoints[i - 1], rotatedPoints[i], ff);
        point.setZorT(z, isDepth);
        return point;
    }

    public StsPoint getPointAtZ(float z, StsSection section, boolean extrapolate)
    {
        StsPoint point = getXYZPointAtZorT(z, extrapolate);
        if (this.onSection == section)
        {
            return point;
        }
        else
        {
            float indexF = section.getLineIndexF(this);
            float vv[] = new float[]{point.getX(), point.getY(), point.getZorT(), indexF};
            return new StsPoint(vv);
        }
    }

    /** perform add action */
    public boolean addLine()
    {
        if (lineVertices.getSize() < 2)
        {
            StsMessageFiles.logMessage("Less than 2 line vertices picked for " +
                getLabel() + ": construction aborted.");
            return false;
        }
        else
        {
            try
            {
                if (onSection != null)
                {
                    addOnSection(onSection);
                }
                currentModel.getSpectrumClass().incrementSpectrumColor("Basic");
            }
            catch (Exception e)
            {
                StsException.outputException("StsLine.addLine() failed.", e, StsException.WARNING);
            }
        }
        currentSectionLine = null;
        return true;
    }

    public boolean addOnSection(StsSection section)
    {
        try
        {
            setOnSection(section);
            section.setLineOnSide(this);
            /*
                     int nVertices = lineVertices.getSize();

                     for (int n = 0; n < nVertices; n++)
                     {
                         StsSurfaceVertex vertex = (StsSurfaceVertex)lineVertices.getElement(n);
                         section.projectVertexToSection(vertex);
                         vertex.setAssociation(this);
                     }
            */
            return projectToSection();
            //            computePoints();
            //            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsLine.addOnSection() failed.",
                e, StsException.WARNING);
            return false;
        }

    }

    public boolean hasConnectedSection()
    {
        if (connectedSections == null)
        {
            return false;
        }
        int nConnectedSections = connectedSections.getSize();
        for (int n = 0; n < nConnectedSections; n++)
        {
            StsSection connectedSection = (StsSection) connectedSections.getElement(n);
            if (connectedSection != null && connectedSection.getIndex() >= 0)
            {
                return true;
            }
        }
        return false;
    }

    public boolean hasSection()
    {
        return onSection != null || (connectedSections != null && connectedSections.getSize() > 0);
    }

    public StsSection getPrevSection()
    {
        if (connectedSections == null)
        {
            return null;
        }
        int nConnectedSections = connectedSections.getSize();
        for (int n = 0; n < nConnectedSections; n++)
        {
            StsSection section = (StsSection) connectedSections.getElement(n);
            if (section.getLastLine() == this)
            {
                return section;
            }
        }
        return null;
    }

    public StsSection getOtherSection(StsSection section)
    {
        if (connectedSections == null)
        {
            return null;
        }
        int nConnectedSections = connectedSections.getSize();
        if (nConnectedSections < 2)
        {
            return null;
        }
        for (int n = 0; n < nConnectedSections; n++)
        {
            StsSection otherSection = (StsSection) connectedSections.getElement(n);
            if (section != otherSection)
            {
                return otherSection;
            }
        }
        return null;
    }

    /** remove a well from the instance list and in the 3d window */
    public boolean delete()
    {
        // don't delete if connected to a section and on a section
        if (hasConnectedSection())
        {
            return false;
        }

        StsObjectRefList.deleteAll(sectionEdgeVertices);
        StsObjectRefList.deleteAll(lineVertices);
        super.delete();
        return true;
    }

    /*
       public void delete(StsSection connectedSection)
       {
           deleteConnectedSection(connectedSection);
           if(!hasConnectedSection())
           {
               if(section != null) section.deleteLine(this); // if line on section, delete from section side
               delete();
           }
       }
    */

    public void deleteSection(StsFracture fracture)
    {
        deleteSection((StsSection) fracture);
    }

    // if well is on this section or connected to this section, delete section from well
    public void deleteSection(StsSection section)
    {
        if (this.onSection == section)
        {
            section.deleteLine(this);
            this.onSection = null;
            sectionSide = NONE;
        }
        else if (connectedSections != null)
        {
            deleteConnectedSection(section);
        }

        if (!hasConnectedSection() && this.onSection != null)
        {
            this.onSection.deleteLine(this); // if line on section, delete from section side
            this.onSection = null;
        }

        delete(); // conditionally delete if no connected sections and not on a section
    }

    /** Remove a surfaceVertex from the appropriate vertices list */
    public void deleteSectionVertex(StsSurfaceVertex vertex)
    {
        if (sectionEdgeVertices == null)
        {
            return;
        }
        sectionEdgeVertices.delete(vertex);
    }

    public void deleteLineVertex(StsSurfaceVertex vertex)
    {
        if (lineVertices == null)
        {
            return;
        }
        lineVertices.delete(vertex);
    }

    public void deleteSurfaceVertices()
    {
        surfaceEdgeVertices = null;
    }

    public boolean computePointsProjectToSection()
    {
        if (!computePoints()) return false;
        return projectToSection();
    }

    public boolean computeExtendedPointsProjectToSection(StsSection section)
    {
        if(section == null) return false;
        return computeExtendedPointsProjectToSection(section.sectionZMin,  section.sectionZMax);
    }

    public boolean computeExtendedPointsProjectToSection(float zMin, float zMax)
    {
        if (!computeExtendedPoints(zMin, zMax)) return false;
        return projectToSection();
    }

    /*
        public boolean computePoints()
        {
            boolean extend = ! (this instanceof StsWell);
            return computeXYZPoints(extend);
        }
    */
    public boolean projectToSection()
    {
        boolean ok = true;

        // projectVerticesToSection();

        if (rotatedPoints == null) return false;
        if (onSection == null) return false;

        if (isVertical && onSection.isVertical())
        {
            StsSectionPoint sectionPoint = new StsSectionPoint(rotatedPoints[0], -1.0f);
            if (onSection.computeNearestPoint(sectionPoint))
            {
                StsPoint nearestPoint = sectionPoint.nearestPoint;
                float colF = sectionPoint.sectionColF;
                for (int n = 0; n < rotatedPoints.length; n++)
                {
                    rotatedPoints[n] = new StsPoint(6, rotatedPoints[n]);
                    rotatedPoints[n].setX(nearestPoint.getX());
                    rotatedPoints[n].setY(nearestPoint.getY());
                    rotatedPoints[n].setZorT(nearestPoint.getZorT());
                    rotatedPoints[n].setF(colF);
                    //                    points[n] = new StsPoint(points[n], colF, isDepth);
                }
            }
            else
            {
                ok = false;
            }
        }
        else
        {
            float guessColF = -1.0f; /** use colF from previous point to speed up search. */
            for (int n = 0; n < rotatedPoints.length; n++)
            {
                StsSectionPoint sectionPoint = new StsSectionPoint(rotatedPoints[n], guessColF);
                if (onSection.computeNearestPoint(sectionPoint))
                {
                    rotatedPoints[n] = new StsPoint(6, rotatedPoints[n]);
                    StsPoint nearestPoint = sectionPoint.nearestPoint;
                    rotatedPoints[n].setX(nearestPoint.getX());
                    rotatedPoints[n].setY(nearestPoint.getY());
                    rotatedPoints[n].setZorT(nearestPoint.getZorT());
                    rotatedPoints[n].setF(sectionPoint.sectionColF);
                    //                   points[n] = new StsPoint(sectionPoint.nearestPoint, sectionPoint.sectionColF, isDepth);
                    guessColF = sectionPoint.sectionColF;
                }
                else
                {
                    ok = false;
                }
            }
            setIsVertical(false);
        }
        return ok;
    }

    public void projectVerticesToSection()
    {
        // project the line vertices to the section
        projectVerticesToSection(lineVertices.getList());
        // using the projected line vertices, spline a new set of points and project each point to the section
        computeExtendedPointsProjectToSection(onSection);
        projectVerticesToSection(sectionEdgeVertices.getList());
        projectVerticesToSection(surfaceEdgeVertices);
    }

    private void projectVerticesToSection(StsObjectList verticesList)
    {
        int nVertices = verticesList.getSize();
        // compute the surfaceVertexes at the intersections of the new splined line and the corresponding surface
        for (int n = 0; n < nVertices; n++)
        {
            StsSurfaceVertex vertex = (StsSurfaceVertex) verticesList.getElement(n);
            onSection.projectVertexToSection(vertex);
            vertex.adjustToGrid();
        }
    }

    /*
       private boolean checkIsVertical()
       {
           int nVertices = (lineVertices == null) ? 0 : lineVertices.getSize();
           if (nVertices < 1) return false;
           StsPoint comparePoint = null;
           for (int i=0; i<nVertices; i++)
           {
               StsSurfaceVertex v = (StsSurfaceVertex)lineVertices.getElement(i);
               if (v==null) continue;
               StsPoint point = v.getPoint();
               if (comparePoint == null) comparePoint = point;
               else if (!comparePoint.sameXY(point)) return false;
           }
           return true;
       }
    */

    public boolean projectRotationAngleChanged()
    {
        return computePoints();
    }

    public boolean computePoints()
    {
        try
        {
            StsPoint[] projectPoints = computeRotatedCoorVertexPoints();
            if (projectPoints.length <= 1) return true;
            StsPoint[] slopes = StsBezier.computeXYZLineSlopes(projectPoints);
            isVertical = checkIsVertical(slopes);
            return computeRotatedPointsFromVertexPoints(projectPoints, slopes);
        }
        catch (Exception e)
        {
            StsException.outputException("Exception in StsLine.computeXYZPoints()", e, StsException.WARNING);
            return false;
        }
    }

    public boolean computeExtendedPoints(float zMin, float zMax)
    {
        try
        {
            StsPoint[] projectPoints = computeRotatedCoorVertexPoints();
            StsPoint[] slopes = StsBezier.computeXYZLineSlopes(projectPoints);
            return computeExtendedXYZPoints(projectPoints, slopes, zMin, zMax);
        }
        catch (Exception e)
        {
            StsException.outputException("Exception in StsLine.computeXYZPoints()", e, StsException.WARNING);
            return false;
        }
    }

    public StsPoint[] computeRotatedCoorVertexPoints()
    {
        try
        {
            if (lineVertices == null)
            {
                return new StsPoint[0];
            }
            int nVertices = lineVertices.getSize();
//            addMeasuredDepth();
            StsProject project = currentModel.getProject();

            StsPoint[] projectPoints = new StsPoint[nVertices];

            // imported line objects like wells and fault sticks are in original unrotated coordinate system:
            // get points in rotated coordinates.  Line objects internally constructed are in rotated coordinates.
            // So if unrotated, compute rotated points.
            if (isVerticesRotated)
            {
                for (int n = 0; n < nVertices; n++)
                {
                    StsSurfaceVertex vertex = (StsSurfaceVertex) lineVertices.getElement(n);
                    projectPoints[n] = new StsPoint(5, vertex.getPoint());
                }
            }
            else
            {
                float dXOrigin = (float) (xOrigin - project.getXOrigin());
                float dYOrigin = (float) (yOrigin - project.getYOrigin());
                for (int n = 0; n < nVertices; n++)
                {
                    StsSurfaceVertex vertex = (StsSurfaceVertex) lineVertices.getElement(n);
                    projectPoints[n] = vertex.getRotatedPoint(project, dXOrigin, dYOrigin);
                }
            }
            return projectPoints;
        }
        catch (Exception e)
        {
            StsException.outputException("StsLine.computeProjectPoints() failed.",
                e, StsException.WARNING);
            return null;
        }
    }

    public void resetSectionLinePoints(StsGridBoundingBox modelBoundingBox)
    {
        if(isVertical) return;
        float zMin = modelBoundingBox.zMin;
        float zInc = modelBoundingBox.zInc;

        int nSectionPoints = modelBoundingBox.getNSlices();
        StsPoint[] newRotatedPoints = new StsPoint[nSectionPoints];
        float z = zMin;
        for(int n = 0; n < nSectionPoints; n++, z += zInc)
            newRotatedPoints[n] = getPointAtZorT(z, true);
        rotatedPoints = newRotatedPoints;
    }

    protected void addMeasuredDepth()
    {
        float measuredDepth = 0.0f;
        StsPoint point0 = getLineVertexPoint(0);
        point0.setM(measuredDepth);
        boolean originalIsDepth = isOriginalDepth();
        float[] xyz1 = point0.getXYZorT(originalIsDepth);
        int nPoints = lineVertices.getSize();
        for (int n = 0; n < nPoints - 1; n++)
        {
            float[] xyz0 = xyz1;
            StsPoint point1 = getLineVertexPoint(n + 1);
            xyz1 = point1.getXYZorT(originalIsDepth);
            measuredDepth += StsMath.distance(xyz0, xyz1, 3);
            point1.setM(measuredDepth);
        }
    }

    protected StsPoint getLineVertexPoint(int index)
    {
        StsSurfaceVertex vertex = (StsSurfaceVertex) lineVertices.getElement(index);
        return vertex.getPoint();
    }

    public StsPoint getRotatedPoint(StsSurfaceVertex vertex)
    {
        StsProject project = currentModel.getProject();
        float dXOrigin = (float) (xOrigin - project.getXOrigin());
        float dYOrigin = (float) (yOrigin - project.getYOrigin());
        return vertex.getRotatedPoint(project, dXOrigin, dYOrigin);
    }

    /**
     * Compute points for this StsLine from vertices.
     * If clipToProject is true, and line extends above and/or below project depth limits,
     * add an interpolated point at project limit(s).
     * If one or both ends do not extend out of project vertically and extend is true,
     * add vertices at project limits as well.
     */
    protected boolean computeRotatedPointsFromVertexPoints(StsPoint[] vertexPoints, StsPoint[] slopes) throws StsException
    {
        float z;
        float f;
        float len;
        int degree = 3;
        StsPoint[] cp = new StsPoint[degree + 1];
        float u, du;

        int nVertexPoints = vertexPoints.length;
        if (nVertexPoints < 1) return false;
        if (nVertexPoints == 1 || slopes == null) return false;

        int nCoors = vertexPoints[0].v.length;

        int nPoints = 0;
        int[] pointsPerInterval = new int[nVertexPoints - 1];
        for (int v = 0; v < nVertexPoints - 1; v++)
        {
            len = vertexPoints[v + 1].getM() - vertexPoints[v].getM();
            int i = (int) (len / maxSplineDistance);
            if (i < 1) i = 1;
            pointsPerInterval[v] = i;
            nPoints += i; // first vertex point plus intermediate points
        }
        nPoints++; // last vertex point

        rotatedPoints = new StsPoint[nPoints];

        // add first point at top of project if extrapolated
        int n = 0;
        for (int v = 0; v < nVertexPoints - 1; v++)
        {
            int ppi = pointsPerInterval[v];
            rotatedPoints[n++] = new StsPoint(vertexPoints[v]);
            if (ppi <= 1) continue;

            cp[0] = vertexPoints[v];
            cp[3] = vertexPoints[v + 1];
            len = cp[3].getM() - cp[0].getM();

            f = len / 3.0f;
            cp[1] = new StsPoint(nCoors);
            cp[1].multByConstantAddPoint(slopes[v], f, vertexPoints[v]);

            f = -len / 3.0f;
            cp[2] = new StsPoint(nCoors);
            cp[2].multByConstantAddPoint(slopes[v + 1], f, vertexPoints[v + 1]);

            u = 0.0f;
            du = 1.0f / ppi;
            for (int i = 0; i < ppi - 1; i++)
            {
                u += du;

                rotatedPoints[n++] = StsBezier.evalBezierCurve(cp, degree, u);

            }

        }
        // add last vertex point
        rotatedPoints[n++] = new StsPoint(vertexPoints[nVertexPoints - 1]);

        return true;
    }

    private boolean checkIsVertical(StsPoint[] slopes)
    {
        for (int n = 0; n < slopes.length; n++)
        {
            if (!StsMath.sameAsTol(slopes[n].getX(), 0.0f, 1.0e-5f)) return false;
            if (!StsMath.sameAsTol(slopes[n].getY(), 0.0f, 1.0e-5f)) return false;
        }
        return true;
    }

    protected boolean computeExtendedXYZPoints(StsPoint[] vertexPoints, StsPoint[] slopes, float zTop, float zBot) throws StsException
    {
        int nVertexTop = 0, nVertexBot = 0;
        float z;
        float f;
        float len;
        int degree = 3;
        StsPoint[] cp = new StsPoint[degree + 1];
        float u, du;

        int nVertexPoints = vertexPoints.length;
        if (nVertexPoints < 1) return false;
        if (nVertexPoints == 1 || slopes == null) return false;

        int nCoors = vertexPoints[0].v.length;

        nVertexTop = getVertexIndexBelow(vertexPoints, zTop);
        if (nVertexTop == nVertexPoints) return false;
        nVertexBot = getVertexIndexAbove(vertexPoints, zBot);
        if (nVertexBot == -1) return false;
        float topVertexZ = vertexPoints[nVertexTop].getZorT();
        float botVertexZ = vertexPoints[nVertexBot].getZorT();

        // put nSplineInterval points between each pair of vertices
        // plus a point above and/or below if extrapolated or interpolated

        int nPoints = 0;
        int nIntervals = nVertexBot - nVertexTop;
        int[] pointsPerInterval = new int[nIntervals];
        for (int v = nVertexTop, interval = 0; v < nVertexBot; v++, interval++)
        {
            len = vertexPoints[v + 1].getM() - vertexPoints[v].getM();
            int i = (int) (len / maxSplineDistance);
            if (i < 1) i = 1;
            pointsPerInterval[interval] = i;
            nPoints += i; // first vertex point plus intermediate points
        }
        nPoints++; // add last vertex point
        if (zTop < topVertexZ)
        {
            nPoints++; // extrapolated top point
        }
        if (zBot > botVertexZ)
        {
            nPoints++; // extrapolated bottom point
        }

        rotatedPoints = new StsPoint[nPoints];

        // add first point at top of project if extrapolated
        int n = 0;
        if (zTop < topVertexZ)
        {
            f = (zTop - topVertexZ) / slopes[nVertexTop].getZorT();
            rotatedPoints[0] = new StsPoint(nCoors);
            StsPoint topVertex = vertexPoints[nVertexTop];
            rotatedPoints[0].multByConstantAddPoint(slopes[nVertexTop], f, topVertex);
            //			points[0].setF(topVertex.getZ()); // top measured depth is vertical distance down to point
            n++;
        }

        for (int v = nVertexTop, interval = 0; v < nVertexBot; v++, interval++)
        {
            int ppi = pointsPerInterval[interval];
            rotatedPoints[n++] = new StsPoint(vertexPoints[v]);
            if (ppi <= 1) continue;

            cp[0] = vertexPoints[v];
            cp[3] = vertexPoints[v + 1];
            len = cp[3].getM() - cp[0].getM();

            f = len / 3.0f;
            cp[1] = new StsPoint(nCoors);
            cp[1].multByConstantAddPoint(slopes[v], f, vertexPoints[v]);

            f = -len / 3.0f;
            cp[2] = new StsPoint(nCoors);
            cp[2].multByConstantAddPoint(slopes[v + 1], f, vertexPoints[v + 1]);

            u = 0.0f;
            du = 1.0f / ppi;
            for (int i = 0; i < ppi - 1; i++)
            {
                u += du;
                rotatedPoints[n++] = StsBezier.evalBezierCurve(cp, degree, u);
            }
        }

        // add last vertex point
        rotatedPoints[n++] = new StsPoint(vertexPoints[nVertexBot]);

        // add extrapolated point at project bottom if extrapolated
        if (zBot > botVertexZ)
        {
            StsPoint botVertex = vertexPoints[nVertexBot];
            f = (zBot - botVertex.getZorT()) / slopes[nVertexBot].getZorT();
            StsPoint botPoint = new StsPoint(nCoors);
            botPoint.multByConstantAddPoint(slopes[nVertexBot], f, botVertex);
            rotatedPoints[n++] = botPoint;
        }
        return true;
    }

    /**
     * Returns interval this z values falls in.
     * Assume the interval above is -1 and interval below is nVertexPoints.
     * So if above or below, return these values.
     */
    private int getVertexIndexBelow(StsPoint[] vertexPoints, float z)
    {
        int nVertexPoints = vertexPoints.length;
        int n = 0;
        for (n = 0; n < nVertexPoints; n++)
            if (vertexPoints[n].getZorT() >= z) return n;
        return nVertexPoints;
    }

    // has not been debugged
    private int getVertexIndexAbove(StsPoint[] vertexPoints, float z)
    {
        int nVertexPoints = vertexPoints.length;
        int n = 0;
        for (n = nVertexPoints - 1; n >= 0; n--)
            if (vertexPoints[n].getZorT() <= z) return n;
        return -1;
    }

    /*
       protected boolean computeXYZPoints(StsPoint[] vertexPoints, StsPoint[] slopes) throws StsException
       {
           StsProject project = currentModel.getProject();
           StsPoint topPoint = null;
           StsPoint botPoint = null;
        float z;
           float zTop, zBot;
           float f;
           float cdz;
           int degree = 3;
           int dim = 3;
           StsPoint[] cp = new StsPoint[degree+1];
           float u, du;

           int nVertexPoints = vertexPoints.length;
           if(nVertexPoints < 1) return false;

           if(nVertexPoints == 1 || slopes == null) return false;

           dZ = project.getZInc();

           float zMin = project.getZMin();
           float botPointZ = vertexPoints[nVertexPoints-1].getZ();
           if (zMin > botPointZ) return false; // all lineVertices are above zMin
           float zMax = project.getZMax();
     float topPointZ = vertexPoints[0].getZ();
           if (zMax < topPointZ) return false; // all lineVertices are below zMax

           int nPointTop = 0;
           for (int n = 0; n < nVertexPoints; n++)
           {
               topPoint = vertexPoints[n];
               topPointZ = topPoint.getZ();
               if (topPointZ >= zMin)
               {
                   nPointTop = n;
                   break;
               }
           }
           if (nPointTop == 0) iTop = project.getIndexAbove(topPointZ);
           else iTop = project.getIndexAbove(Math.min(zMin,topPointZ));
           zTop = iTop*dZ;

           int nPointBottom = nVertexPoints - 1;
           for (int n = nPointTop; n < nVertexPoints; n++)
           {
               botPoint = vertexPoints[n];
               botPointZ = botPoint.getZ();
               if (botPointZ == zMax)
               {
                   nPointBottom = n;
                   break;
               }
               if (botPointZ > zMax)
               {
                   if (n == nPointTop)
                   {
                       nPointBottom = n;
                       break;
                   }
                   nPointBottom = n - 1;
                   botPoint = vertexPoints[n-1];
                   botPointZ = botPoint.getZ();
                   break;
               }
           }
           if (nPointTop == nPointBottom)
           {
               iBot = project.getIndexBelow(Math.max(zMax, botPointZ));
           }
           else
           {
               iBot = project.getIndexBelow(Math.min(zMax, botPointZ));
           }
           zBot = iBot*dZ;

           int nPoints = iBot - iTop + 1;
           points = new StsPoint[nPoints];

           z = zTop;
           int nStart = 0;
           if (zTop == topPointZ)
           {
               points[0] = topPoint;
               nStart++;
               z += dZ;
           }
           else if (zTop<topPointZ)
           {
               if (nPointTop==0)  // first point is extrapolated
               {
                   f = (zTop - topPoint.getZ())/slopes[nPointTop].getZ();
                   points[0] = new StsPoint();
                   points[0].multByConstantAddPoint(slopes[nPointTop], f, topPoint);
                   nStart++;
                   z += dZ;
               }
               else nPointTop--;
           }

           int nEnd = nPoints-1;
           if (zBot == botPointZ)
           {
               points[nEnd] = botPoint;
               nEnd--;
           }
           else if (zBot > botPointZ)
           {
               if (nPointBottom==nPoints-1) // last point is extrapolated
               {
                   f = (zBot - botPoint.getZ())/slopes[nPointBottom].getZ();
                   points[nEnd] = new StsPoint();
                   points[nEnd].multByConstantAddPoint(slopes[nPointBottom],
                           f, botPoint);
                   nEnd--;
               }
               else nPointBottom++;
           }

           int v = nPointTop;
           u = 0.0f;
           du = 1.0f;
           float nextVertexZ = vertexPoints[v+1].getZ();

           for (int n = nStart; n <= nEnd; n++)
           {
               if(n == nStart || z > nextVertexZ)
               {
                   while (v < nPointBottom-2 && z > nextVertexZ)
                   {
                       v++;
                       nextVertexZ = vertexPoints[v+1].getZ();
                   }
                   cp[0] = vertexPoints[v];
                   cp[3] = vertexPoints[v+1];

                   cdz = cp[3].getZ() - cp[0].getZ();

                   f = cdz/(3.0f*slopes[v].getZ());
                   cp[1] = new StsPoint();
                   cp[1].multByConstantAddPoint(slopes[v], f, vertexPoints[v]);

                   f = -cdz/(3.0f*slopes[v+1].v[2]);
                   cp[2] = new StsPoint();
                   cp[2].multByConstantAddPoint(slopes[v+1], f, vertexPoints[v+1]);

                   u = (z - cp[0].getZ())/cdz;
                   du = dZ/cdz;
               }

               points[n] = StsBezier.evalBezierCurve(cp, degree, u);
//			points[n].print("point: ", n);
               u += du;
               z += dZ;
           }
//        points[n].print("point: ", nPoints-1);
           return true;
    }
    */
    /**
     * Try an iterative top-down method to find well-grid intersection.  If this fails,
     * go to a methodical top-down search.
     */

    public StsGridPoint computeGridIntersect(StsXYSurfaceGridable grid)
    {
        if (rotatedPoints == null)
        {
            return null;
        }
        if (isVertical || rotatedPoints.length < 2)
        {
            return computeVerticalGridIntersect(grid);
        }
        else
        {
            return computeGridIntersect(grid, grid.getZMin());
        }
        //            return computeGridIntersect(grid, grid.getZMin(), grid.getZMax());
    }

    // from this startZ iterate to final intersection
    public StsGridPoint computeGridIntersect(StsXYSurfaceGridable grid, float startZ)
    {
        int indexAbove, indexBelow;
        boolean aboveOK, belowOK;
        StsGridPoint point;
        try
        {
            int nPoints = rotatedPoints.length;

            if (startZ < rotatedPoints[0].getZorT())
            {
                indexAbove = -1;
            }
            else if (startZ > rotatedPoints[nPoints - 1].getZorT())
            {
                indexAbove = nPoints - 1;
            }
            else
            {
                indexAbove = 0;
                float z1 = rotatedPoints[0].getZorT();
                for (int n = 1; n < nPoints; n++)
                {
                    float z0 = z1;
                    z1 = rotatedPoints[n].getZorT();
                    if (startZ >= z0 && startZ <= z1)
                    {
                        indexAbove = n - 1;
                        break;
                    }
                }
            }
            indexBelow = indexAbove + 1;
            aboveOK = indexAbove >= 0;
            belowOK = indexBelow < nPoints;

            while (aboveOK || belowOK)
            {
                if (aboveOK)
                {
                    point = computeGridIntervalIntersect(grid, indexAbove + 1);
                    if (point != null)
                    {
                        return point;
                    }
                    indexAbove--;
                    if (indexAbove < 0)
                    {
                        aboveOK = false;
                    }
                }
                if (belowOK)
                {
                    point = computeGridIntervalIntersect(grid, indexBelow);
                    if (point != null)
                    {
                        return point;
                    }
                    indexBelow++;
                    if (indexBelow >= nPoints)
                    {
                        belowOK = false;
                    }
                }
            }

            // didn't intersect: use top or bot point, whichever is closer
            /*
                        StsPoint topPoint = points[0];
                        float gridZ = grid.interpolateBilinearZ(topPoint, true, false);
                        float topDZ = Math.abs(gridZ - topPoint.getZorT());
                        StsPoint botPoint = points[nPoints - 1];
                        gridZ = grid.interpolateBilinearZ(botPoint, true, false);
                        float botDZ = Math.abs(gridZ - botPoint.getZorT());
                        if (topDZ <= botDZ)
                        {
                            indexBelow = 0;
                        }
                        else
                        {
                            return botPoint;
                        }
            */
            // didn't find intersection along wellLine: try extrapolations above and below

            point = computeGridIntervalIntersect(grid, -1);
            if (point != null) return point;
            point = computeGridIntervalIntersect(grid, nPoints);
            if (point != null) return point;
            StsException.systemError("StsLine.computeGridIntersect() failed." +
                " Fault: " + getLabel() + " Grid: " + grid.getLabel());
            return null;
            
        }
        catch (Exception e)
        {
            StsException.outputException("StsLine.computeGridIntersect() failed." +
                " Fault: " + getLabel() + " Grid: " + grid.getLabel(),
                e, StsException.WARNING);
            return null;
        }
    }

    private StsGridPoint computeGridIntervalIntersect(StsXYSurfaceGridable grid, int indexBelow)
    {
        StsPoint point0, point1;

        if (rotatedPoints == null)
        {
            return null;
        }

        if (indexBelow <= 0)
        {
            point0 = getXYZPointAtZorT(grid.getZMin(), true);
            point1 = rotatedPoints[0].getXYZorTPoint();
        }
        else if (indexBelow > rotatedPoints.length - 1)
        {
            point0 = rotatedPoints[rotatedPoints.length - 1].getXYZorTPoint();
            point1 = getXYZPointAtZorT(grid.getZMax(), true);
        }
        else
        {
            point0 = rotatedPoints[indexBelow - 1].getXYZorTPoint();
            point1 = rotatedPoints[indexBelow].getXYZorTPoint();
        }

        StsGridCrossings gridCrossings = new StsGridCrossings(point0, point1, grid);
        return gridCrossings.getGridIntersection(grid);
    }

    public StsSurfaceVertex getSectionEdgeVertex(StsXYSurfaceGridable surface)
    {
        StsSurfaceVertex vertex = null;
        try
        {
            if(sectionEdgeVertices != null)
                vertex = getEdgeVertex(sectionEdgeVertices.getList(), surface);
            if (vertex != null) return vertex;
            vertex = computeSectionEdgeSurfaceVertex(surface);
            insertSectionEdgeVertex(vertex);
            return vertex;
        }
        catch (Exception e)
        {
            StsException.outputException("StsLine.getSurfaceVertex() failed.",
                e, StsException.WARNING);
            return null;
        }
    }

    public StsSurfaceVertex getSurfaceEdgeVertex(StsXYSurfaceGridable surface)
    {
        StsSurfaceVertex vertex = null;
        try
        {
            if(surfaceEdgeVertices != null)
                vertex = getEdgeVertex(surfaceEdgeVertices, surface);
            if (vertex != null) return vertex;
            vertex = computeSurfaceEdgeVertex(surface);
            insertSurfaceEdgeVertex(vertex);
            return vertex;
        }
        catch (Exception e)
        {
            StsException.outputException("StsLine.getSurfaceVertex() failed.",
                e, StsException.WARNING);
            return null;
        }
    }

    public StsSurfaceVertex getEdgeVertex(StsObjectList edgeVertices, StsXYSurfaceGridable surface)
    {
        if (edgeVertices == null) return null;
        int nVertices = edgeVertices.getSize();
        for (int n = 0; n < nVertices; n++)
        {
            StsSurfaceVertex vertex = (StsSurfaceVertex) edgeVertices.getElement(n);
            if (vertex.getSurface() == surface)
            {
                return vertex;
            }
        }
        return null;
    }

    public boolean hasSectionEdgeVertex(StsSurface surface)
    {
        if (sectionEdgeVertices == null)
        {
            return false;
        }
        int nVertices = sectionEdgeVertices.getSize();
        for (int n = 0; n < nVertices; n++)
        {
            StsSurfaceVertex vertex = (StsSurfaceVertex) sectionEdgeVertices.getElement(n);
            if (vertex.getSurface() == surface)
            {
                return true;
            }
        }
        return false;
    }

    public StsSurfaceVertex computeSectionEdgeSurfaceVertex(StsXYSurfaceGridable surface)
    {
        return computeSurfaceVertex(surface, true);
    }

    public StsSurfaceVertex computeSurfaceEdgeVertex(StsXYSurfaceGridable surface)
    {
        return computeSurfaceVertex(surface, false);
    }

    public StsSurfaceVertex computeSurfaceVertex(StsXYSurfaceGridable surface, boolean persistent)
    {
        StsGridPoint gridPoint = computeGridIntersect(surface);
        if (gridPoint == null)
        {
            return null;
        }
        surface.interpolateBilinearZ(gridPoint, true, true);
        return new StsSurfaceVertex(gridPoint.getPoint(), this, surface, null, persistent);
    }
/*
    public void adjustSurfaceVertex(StsXYSurfaceGridable surface)
    {
        StsGridPoint gridPoint = computeGridIntersect(surface);
        if (gridPoint == null) return;
        surface.interpolateBilinearZ(gridPoint, true, true);
        StsSurfaceVertex surfaceVertex = getSurfaceVertex(surface);
        surfaceVertex.setPoint(gridPoint.getPoint());
    }
*/
    // use this method in model construction as it has only x, y, and z
/*
    private StsSurfaceVertex computeModelSurfaceVertex(StsSurface surface, StsBlock block)
    {
        StsGridPoint gridPoint = computeGridIntersect(surface);
        if (gridPoint == null)
        {
            return null;
        }
        surface.interpolateBilinearZ(gridPoint, true, true);
        //        StsPoint point = new StsPoint(gridPoint);
        StsPoint point = new StsPoint(gridPoint.getXYZorT(isDepth));
        //        getUnrotatedRelativeXYFromRotatedXY(point);
        return new StsSurfaceVertex(point, this, surface, block, true);
    }
*/
    /*
       public StsSurfaceVertex addLineVertex(StsPoint point)
       {
           return addLineVertex(point, true);
       }
    */
    /** Add a vertex to the list of vertices in increasing order of z. */
    public StsSurfaceVertex addLineVertex(StsPoint point, boolean computePath, boolean extend)
    {
        StsSurfaceVertex vertex;

        if (point == null)
        {
            StsMessageFiles.logMessage("Can't add a null point.");
            return null;
        }
        float z = point.getZorT();
        if (point.getLength() < 5)
        {
            if (isDepth)
            {
                point = new StsPoint(point.getX(), point.getY(), z, 0.0f, 0.0f);
            }
            else
            {
                point = new StsPoint(point.getX(), point.getY(), 0.0f, 0.0f, z);
            }
        }
        StsProject project = currentModel.getProject();
        if (z < project.getZorTMin())
        {
            StsMessageFiles.logMessage("Can't add z: " + z + ". Less than project min: "
                + project.getZorTMin());
            return null;
        }
        else if (z > project.getZorTMax())
        {
            StsMessageFiles.logMessage("Can't add z: " + z + ". Greater than project max: "
                + project.getZorTMax());
            return null;
        }

        /** Check if picked point seems like a large lateral pick (mispick or start
         *  of next well. Query the user if he wants it.
         */
        if (lineVertices.getSize() > 0)
        {
            StsSurfaceVertex lastVertex = (StsSurfaceVertex) lineVertices.getLast();
            StsPoint lastPoint = lastVertex.getPoint();
            float dH = point.distanceSquaredType(StsPoint.XY, lastPoint);
            float dV = point.distanceSquaredType(StsPoint.Z, lastPoint);

            if (checkExcursion && dV > 0.0f && dH > 2.0f * dV)
            {
                if (!StsYesNoDialog.questionValue(currentModel.win3d, "Lateral excursion of this pick seems large.\n" +
                    "Do you wish to include this point?"))
                {
                    return null;
                }
            }
        }

        /** If we don't have a color yet, get it now */
        if (stsColor == null)
        {
            stsColor = new StsColor(currentModel.getSpectrumClass().getCurrentSpectrumColor("Basic"));
        }

        vertex = new StsSurfaceVertex(point, this);
        // measured depth is used as dimensionaless parameter in splining curve shape, so set Z or T as value
        point.setM(z);
        //        float dZMax = 0.0f;
        //        try { dZMax = 2.0f * project.getZInc(); }
        //        catch (Exception e) { return null; }
        if (!insertLineVertex(vertex, computePath, 1.0f, false)) // use tolerance
        {
            lineVertices.delete(vertex);
            vertex = null; // garbage collect
            return null;
        }
        return vertex;
    }

    public boolean insertLineVertex(StsSurfaceVertex vertex)
    {
        return insertLineVertex(vertex, false, StsParameters.smallFloat, false);
    }

    public boolean insertLineVertex(StsSurfaceVertex vertex, boolean computePath)
    {
        return insertLineVertex(vertex, computePath, StsParameters.smallFloat, false);
    }

    public boolean insertLineVertex(StsSurfaceVertex vertex, boolean computePath,
                                    float minZTolerance)
    {
        return insertLineVertex(vertex, computePath, minZTolerance, false);
    }

    public boolean insertLineVertex(StsSurfaceVertex vertex, boolean computePath, float minZTolerance, boolean replaceWithinTolerance)
    {
        if (vertex == null)
        {
            return false;
        }
        if (lineVertices.contains(vertex))
        {
            return false;
        }

        StsPoint p = vertex.getPoint();
        if (p == null)
        {
            return false;
        }
        StsProject proj = currentModel.getProject();
        if (proj == null)
        {
            return false;
        }

        float z = p.getZorT();
        if (z < proj.getZorTMin() || z > proj.getZorTMax())
        {
            return false;
        }

        float lastZ = -StsParameters.largeFloat;
        if (z == lastZ)
        {
            return false;
        }

        vertex.setAssociation(this);

        if (onlyMonotonic())
        {
            StsSurfaceVertex lineVertex = (StsSurfaceVertex) lineVertices.getFirst();
            while (lineVertex != null)
            {
                float nextZ = lineVertex.getPoint().getZorT();
                if (z == lastZ)
                {
                    return false;
                }

                if (z >= lastZ && z <= nextZ)
                {
                    if (z - lastZ < minZTolerance || nextZ - z < minZTolerance)
                    {
                        if (replaceWithinTolerance)
                        //                    if (replaceWithinTolerance && !lineVertex.usedInLine(currentModel))
                        {
                            try
                            {
                                lineVertices.delete(lineVertex);
                                float vertexZ = vertex.getPoint().getZorT();
                                int nVertices = lineVertices.getSize();
                                for (int i = 0; i < nVertices; i++)
                                {
                                    StsSurfaceVertex v = (StsSurfaceVertex) lineVertices.getElement(i);
                                    if (v.getPoint().getZorT() > vertexZ)
                                    {
                                        lineVertices.insertBefore(v, vertex);
                                        if (computePath)
                                        {
                                            computePointsProjectToSection();
                                        }
                                        return true;
                                    }
                                }
                                lineVertices.add(vertex);
                                if (computePath)
                                {
                                    computePointsProjectToSection();
                                }
                                return true;
                            }
                            catch (Exception e)
                            {
                                return false;
                            }
                        }
                        return false;
                    }
                    lineVertices.insertBefore(lineVertex, vertex);
                    if (computePath)
                    {
                        computePointsProjectToSection();
                    }
                    return true;
                }
                lastZ = nextZ;
                lineVertex = (StsSurfaceVertex) lineVertices.getNext();
            }
        }
        lineVertices.add(vertex);
        if (computePath)
        {
            computePointsProjectToSection();
        }
        return true;
    }

    public boolean onlyMonotonic() { return true; }

    public boolean insertSectionEdgeVertex(StsSurfaceVertex vertex)
    {
        if(sectionEdgeVertices == null)
            sectionEdgeVertices = StsObjectRefList.constructor(2, 2, "sectionEdgeVertices", this);
        return insertVertex(vertex, sectionEdgeVertices.getList());
    }

    public boolean insertSurfaceEdgeVertex(StsSurfaceVertex vertex)
    {
        if(surfaceEdgeVertices == null)
            surfaceEdgeVertices = new StsObjectList(2, 2);
        return insertVertex(vertex, surfaceEdgeVertices);
    }

    public boolean insertVertex(StsSurfaceVertex vertex, StsObjectList surfaceVertices)
    {
        if (vertex == null) return false;
        if (surfaceVertices.contains(vertex))
            return false;

        StsPoint p = vertex.getPoint();
        if (p == null)
            return false;
        StsProject proj = currentModel.getProject();
        if (proj == null)
            return false;

        float z = p.getZorT();
        if (z < proj.getZorTMin() || z > proj.getZorTMax())
            return false;

        float lastZ = -StsParameters.largeFloat;
        if (z == lastZ)
            return false;

        vertex.setAssociation(this);

        StsSurfaceVertex lastVertex, nextVertex;
        float nextZ;

        int nSurfaceVertices = surfaceVertices.getSize();
        if (nSurfaceVertices == 0)
            surfaceVertices.add(vertex);
        else if (nSurfaceVertices == 1)
        {
            nextVertex = (StsSurfaceVertex) surfaceVertices.getElement(0);
            nextZ = nextVertex.getPoint().getZorT();
            if (z <= nextZ)
                surfaceVertices.insertBefore(0, vertex);
            else
                surfaceVertices.add(vertex);
        }
        else
        {
            nextVertex = (StsSurfaceVertex) surfaceVertices.getElement(0);
            nextZ = nextVertex.getPoint().getZorT();

            for (int n = 1; n < nSurfaceVertices; n++)
            {
                lastZ = nextZ;
                nextVertex = (StsSurfaceVertex) surfaceVertices.getElement(n);
                nextZ = nextVertex.getPoint().getZorT();
                if (z >= lastZ && z <= nextZ)
                {
                    surfaceVertices.insertBefore(n, vertex);
                    return true;
                }
            }
            surfaceVertices.add(vertex);
        }
        return true;
    }

    public StsGridPoint computeVerticalGridIntersect(StsXYSurfaceGridable grid)
    {
        StsPoint point = rotatedPoints[0].copy();
        StsGridPoint gridPoint = new StsGridPoint(point, grid);
        float z = grid.interpolateBilinearZ(gridPoint, true, true);
        if (z == nullValue)
        {
            StsException.systemError("StsLine.computeVerticalGridIntersect() failed." +
                " For well: " + getLabel());
            return null;
        }
        return gridPoint;
    }

    // 3D display routines

    /** display this well */

    static public void displayClass(StsGLPanel3d glPanel3d, StsClass instanceList)
    {
        int nInstances = instanceList.getSize();
        for (int n = 0; n < nInstances; n++)
        {
            StsLine line = (StsLine) instanceList.getElement(n);
            line.display(glPanel3d, false);
        }
    }

    public void display(StsGLPanel3d glPanel3d)
    {
        display(glPanel3d, false);
    }

    public void display(StsGLPanel3d glPanel3d, boolean displayName)
    {
        if (glPanel3d == null)
        {
            return;
        }
        if (isVisible)
        {
            display(glPanel3d, highlighted, getName(), rotatedPoints);
        }
    }

    public void display(StsGLPanel3d glPanel3d, String name, StsPoint[] points)
    {
        if (glPanel3d == null)
        {
            return;
        }
        if (isVisible)
        {
            display(glPanel3d, highlighted, name, points);
        }
    }

    public void display(StsGLPanel3d glPanel3d, boolean highlighted, String name, StsPoint[] points)
    {
        display(glPanel3d, highlighted, name, points, false);
    }

    public void display(StsGLPanel3d glPanel3d, boolean highlighted, String name, StsPoint[] points, boolean drawDotted)
    {
        if (glPanel3d == null) return;
        if (!isVisible) return;
        GL gl = glPanel3d.getGL();
        if (gl == null) return;

        if (points == null)
        {
            if (drawVertices) displayVertices(glPanel3d);
            return;
        }

        if(!currentModel.getProject().supportsZDomain(zDomainSupported)) return;

        int nPoints = points.length;
        int zIndex = getPointsZIndex();
        StsPoint topPoint = points[0];
        StsPoint botPoint = points[nPoints - 1];
        StsCropVolume cropVolume = currentModel.getProject().cropVolume;
        boolean isZCropped = cropVolume.isZCropped();
        if (isZCropped)
        {
            float cropMinZ = cropVolume.zMin;
            int min = 0;
            for (int n = 0; n < nPoints; n++)
            {
                if (points[n].getZorT() >= cropMinZ)
                {
                    min = n;
                    topPoint = getLinePointAtZ(cropMinZ, false);
                    break;
                }
            }
            float cropMaxZ = cropVolume.zMax;
            int max = nPoints - 1;
            for (int n = nPoints - 1; n >= 0; n--)
            {
                if (points[n].getZorT() <= cropMaxZ)
                {
                    max = n;
                    botPoint = getLinePointAtZ(cropMaxZ, false);
                    break;
                }
            }
            drawLine(gl, stsColor, highlighted, points, min, max, topPoint, botPoint, zIndex, drawDotted);
        }
        else
        {
            drawLine(gl, stsColor, highlighted, points, zIndex, drawDotted);
        }

        if (name != null)
        {
            displayName(glPanel3d, name, topPoint, botPoint, zIndex);
        }
        if (drawVertices)
        {
            displayVertices(glPanel3d);
        }
        return;
    }

    private void drawLine(GL gl, StsColor stsColor, boolean highlighted, StsPoint[] points, int zIndex, boolean drawDotted)
    {
        int nPoints = points.length;
        drawLine(gl, stsColor, highlighted, points, 0, nPoints - 1, null, null, zIndex, drawDotted);
    }

    public void drawLine(GL gl, StsColor stsColor, boolean highlighted, StsPoint[] points, int min, int max, StsPoint topPoint, StsPoint botPoint, int zIndex, boolean drawDotted)
    {
        if (drawDotted)
            StsGLDraw.drawDottedLine(gl, stsColor, highlighted, points, min, max, topPoint, botPoint, zIndex);
        else
            StsGLDraw.drawLine(gl, stsColor, highlighted, points, min, max, topPoint, botPoint, zIndex);
    }

    public int getPointsZIndex()
    {
        return currentModel.getProject().getPointsZIndex();
    }

    public void displayName(StsGLPanel3d glPanel3d, String name, StsPoint topPoint, StsPoint botPoint, int zIndex)
    {
        float x, y, z;

        GL gl = glPanel3d.getGL();
        if (gl == null)
        {
            return;
        }

        gl.glDisable(GL.GL_LIGHTING);
        stsColor.setGLColor(gl);
        x = topPoint.getX();
        y = topPoint.getY();
        z = topPoint.getZorT();
        StsGLDraw.fontHelvetica12(gl, x, y, z, name);
        x = botPoint.getX();
        y = botPoint.getY();
        z = botPoint.getZorT();
        StsGLDraw.fontHelvetica12(gl, x, y, z, name);
        gl.glEnable(GL.GL_LIGHTING);
    }

    private void displayName(StsGLPanel3d glPanel3d, String name, int min, int max)
    {
        if (rotatedPoints == null)
        {
            return;
        }
        GL gl = glPanel3d.getGL();
        if (gl == null)
        {
            return;
        }

        gl.glDisable(GL.GL_LIGHTING);
        stsColor.setGLColor(gl);
        StsGLDraw.fontHelvetica12(gl, rotatedPoints[min].v, name);
        StsGLDraw.fontHelvetica12(gl, rotatedPoints[max].v, name);
        gl.glEnable(GL.GL_LIGHTING);
    }

    public void display2d(StsGLPanel3d glPanel3d, boolean displayName, int dirNo,
                          float dirCoordinate, boolean axesFlipped, boolean xAxisReversed, boolean yAxisReversed)
    {
        if (glPanel3d == null)
        {
            return;
        }
        if (!isVisible)
        {
            return;
        }

        if (rotatedPoints == null) return;
        if (rotatedPoints.length < 2) return;
        GL gl = glPanel3d.getGL();
        if (gl == null)
        {
            return;
        }

        display2d(glPanel3d, rotatedPoints.length, displayName, dirNo, dirCoordinate, axesFlipped, xAxisReversed, yAxisReversed);
    }

    public void display2d(StsGLPanel3d glPanel3d, int displayToIndex, boolean displayName, int dirNo,
                          float dirCoordinate, boolean axesFlipped, boolean xAxisReversed, boolean yAxisReversed)
    {
        GL gl = glPanel3d.getGL();
        stsColor.setGLColor(gl);
        if (displayName) displayName2d(gl, getName(), dirNo, axesFlipped);
        gl.glDisable(GL.GL_LIGHTING);
        glPanel3d.setViewShift(gl, StsGraphicParameters.edgeShift);
        displayLine2d(gl, stsColor, displayToIndex, dirNo, dirCoordinate, axesFlipped, xAxisReversed, yAxisReversed);
        glPanel3d.resetViewShift(gl);
        gl.glEnable(GL.GL_LIGHTING);
    }

    private void displayLine2d(GL gl, StsColor color, int toIndex, int dirNo, float dirCoordinate, boolean axesFlipped, boolean xAxisReversed, boolean yAxisReversed)
    {
        if (highlighted)
        {
            gl.glLineWidth(StsGraphicParameters.well3dLineWidthHighlighted);
        }
        else
        {
            gl.glLineWidth(StsGraphicParameters.well3dLineWidth);
        }

        gl.glLineStipple(1, StsGraphicParameters.dottedLine);

        int verticalIndex = StsPoint.getVerticalIndex();
        switch (dirNo)
        {
            case 0:
                displayLine2d(gl, color, rotatedPoints, toIndex, 1, verticalIndex, 0, yAxisReversed, dirCoordinate);
                break;
            case 1:
                displayLine2d(gl, color, rotatedPoints, toIndex, 0, verticalIndex, 1, xAxisReversed, dirCoordinate);
                break;
            case 2:
                if (!axesFlipped)
                {
                    displayLine2d(gl, color, rotatedPoints, toIndex, 0, 1, verticalIndex, false, dirCoordinate);
                }
                else
                {
                    displayLine2d(gl, color, rotatedPoints, toIndex, 1, 0, verticalIndex, false, dirCoordinate);
                }
        }
    }

    static final int FRONT = 1;
    static final int BACK = -1;
    static final int END = 0;

    private void displayLine2d(GL gl, StsColor color, StsPoint[] points, int toIndex, int nXAxis2d, int nYAxis2d, int nDepthAxis, boolean axisReversed, float cursorDepth)
    {
        int[] range = new int[]{-1, -1};
        while (true)
        {
            int inFront = isPointsInFrontOrBack(points, nXAxis2d, nYAxis2d, nDepthAxis, axisReversed, range, cursorDepth);
            if(toIndex < range[1])
                range[1] = toIndex;
            switch (inFront)
            {
                case FRONT:
                    StsGLDraw.drawLine2d(gl, color, false, points, range[0], range[1], nXAxis2d, nYAxis2d);
                    break;
                case BACK:
                    StsGLDraw.drawDottedLine2d(gl, color, false, points, range[0], range[1], nXAxis2d, nYAxis2d);
                    break;
                default:
                    return;
            }

        }
    }

    /**
     * draw 2d line between last XY and XY. Line is solid if in front of cursor, dotted if behind.
     * Line is solid if depth is greater than cursorDepth unless axis is reversed; dotted otherwise.
     *
     * @param cursorDepth  depth of cursor
     * @param axisReversed depth direction is reversed
     */

    private int isPointsInFrontOrBack(StsPoint[] rotatedPoints, int nXAxis2d, int nYAxis2d, int nDepthAxis, boolean axisReversed, int[] range, float cursorDepth)
    {
        int nPoints = rotatedPoints.length;
        int min = range[0] + 1;
        if (min >= nPoints) return END;

        int lastInFront = FRONT, inFront;
        float depth = rotatedPoints[min].v[nDepthAxis];
        inFront = getInFront(axisReversed, depth, cursorDepth);
        int max = min + 1;
        for (; max < nPoints; max++)
        {
            lastInFront = inFront;
            inFront = getInFront(axisReversed, depth, cursorDepth);
            if (lastInFront != inFront)
            {
                range[0] = min;
                range[1] = max;
                return lastInFront;
            }
        }
        range[0] = min;
        range[1] = nPoints - 1;
        return lastInFront;
    }

    private int getInFront(boolean axisReversed, float depth, float cursorDepth)
    {
        boolean inFront = axisReversed ? depth >= cursorDepth : depth <= cursorDepth;
        if (inFront)
            return FRONT;
        else
            return BACK;
    }

    private void drawLinePoint2d(GL gl, float lastX, float lastY, float lastDepth, float x, float y, float depth,
                                 float cursorDepth, boolean axisReversed)
    {
        boolean lastInFront, inFront;

        if (axisReversed)
        {
            lastInFront = lastDepth <= cursorDepth;
            inFront = depth <= cursorDepth;
        }
        else
        {
            lastInFront = lastDepth >= cursorDepth;
            inFront = depth >= cursorDepth;
        }
        if (lastInFront == inFront || depth == lastDepth)
        {
            if (!inFront)
            {
                gl.glEnable(GL.GL_LINE_STIPPLE);
            }
            gl.glBegin(GL.GL_LINES);
            gl.glVertex2f(lastX, lastY);
            gl.glVertex2f(x, y);
            gl.glEnd();
            if (!inFront)
            {
                gl.glDisable(GL.GL_LINE_STIPPLE);
            }
        }
        else
        {
            float f = (cursorDepth - lastDepth) / (depth - lastDepth);
            float cursorX = lastX + f * (x - lastX);
            float cursorY = lastY + f * (y - lastY);
            if (inFront)
            {
                gl.glBegin(GL.GL_LINES);
                gl.glVertex2f(lastX, lastY);
                gl.glVertex2f(cursorX, cursorY);
                gl.glEnd();

                gl.glEnable(GL.GL_LINE_STIPPLE);
                gl.glBegin(GL.GL_LINES);
                gl.glVertex2f(cursorX, cursorY);
                gl.glVertex2f(x, y);
                gl.glEnd();
                gl.glDisable(GL.GL_LINE_STIPPLE);
            }
            else
            {
                gl.glEnable(GL.GL_LINE_STIPPLE);
                gl.glBegin(GL.GL_LINES);
                gl.glVertex2f(lastX, lastY);
                gl.glVertex2f(cursorX, cursorY);
                gl.glEnd();
                gl.glDisable(GL.GL_LINE_STIPPLE);

                gl.glBegin(GL.GL_LINES);
                gl.glVertex2f(cursorX, cursorY);
                gl.glVertex2f(x, y);
                gl.glEnd();
            }
            StsGLDraw.drawPoint2d(x, y, gl, 4);
        }
    }

    private void displayName2d(GL gl, String name, int dirNo, boolean axesFlipped)
    {
        float[] xyzTop = rotatedPoints[0].v;
        float[] xyzBot = rotatedPoints[rotatedPoints.length - 1].v;

        switch (dirNo)
        {
            case 0:
                StsGLDraw.fontHelvetica12(gl, xyzTop[1], xyzTop[2], name);
                StsGLDraw.fontHelvetica12(gl, xyzBot[1], xyzBot[2], name);
                break;
            case 1:
                StsGLDraw.fontHelvetica12(gl, xyzTop[0], xyzTop[2], name);
                StsGLDraw.fontHelvetica12(gl, xyzBot[0], xyzBot[2], name);
                break;
            case 2:
                if (!axesFlipped)
                {
                    StsGLDraw.fontHelvetica12(gl, xyzTop[0], xyzTop[1], name);
                    StsGLDraw.fontHelvetica12(gl, xyzBot[0], xyzBot[1], name);
                }
                else
                {
                    StsGLDraw.fontHelvetica12(gl, xyzTop[1], xyzTop[0], name);
                    StsGLDraw.fontHelvetica12(gl, xyzBot[1], xyzBot[0], name);
                }
                break;
        }
    }

    public void displayLabel2d(GL gl, float[] xyz, String name, int dirNo, boolean axesFlipped)
    {
        switch (dirNo)
        {
            case 0:
                StsGLDraw.fontHelvetica12(gl, xyz[1], xyz[2], name);
                StsGLDraw.fontHelvetica12(gl, xyz[1], xyz[2], name);
                break;
            case 1:
                StsGLDraw.fontHelvetica12(gl, xyz[0], xyz[2], name);
                StsGLDraw.fontHelvetica12(gl, xyz[0], xyz[2], name);
                break;
            case 2:
                if (!axesFlipped)
                {
                    StsGLDraw.fontHelvetica12(gl, xyz[0], xyz[1], name);
                    StsGLDraw.fontHelvetica12(gl, xyz[0], xyz[1], name);
                }
                else
                {
                    StsGLDraw.fontHelvetica12(gl, xyz[1], xyz[0], name);
                    StsGLDraw.fontHelvetica12(gl, xyz[1], xyz[0], name);
                }
                break;
        }
    }
    public void displayVertices(StsGLPanel3d glPanel3d)
    {
        //if (win3d == null) return;
        int nVertices = (lineVertices == null) ? 0 : lineVertices.getSize();
        StsPoint[] rotatedVertexPoints = computeRotatedCoorVertexPoints();
        for (int i = 0; i < nVertices; i++)
        {
            StsSurfaceVertex v = (StsSurfaceVertex) lineVertices.getElement(i);
            float[] xyz = rotatedVertexPoints[i].getXYZorT();
            if (v == selectedVertex)
            {
                StsGLDraw.drawPoint(xyz, StsColor.WHITE, glPanel3d, StsGraphicParameters.vertexDotWidthHighlighted);
            }
            else
            {
                StsGLDraw.drawPoint(xyz, StsColor.WHITE, glPanel3d, StsGraphicParameters.vertexDotWidth);
            }
        }
    }

    /*
         public boolean intersectsCursor(StsGLPanel3d glPanel3d, int dirNo)
         {
             StsBoundingBox wellBoundingBox = getBoundingBox();
             return glPanel3d.getCursor3d().cursorIntersected(wellBoundingBox, dirNo);
         }
    */
    public void pick(GL gl, StsGLPanel glPanel)
    {
        if(!currentModel.getProject().canDisplayZDomain(zDomainSupported)) return;
        int zIndex = getPointsZIndex();
        if (rotatedPoints == null)
        {
            return;
        }
        if(!currentModel.getProject().supportsZDomain(zDomainSupported))
        {
            return;
        }
        StsGLDraw.pickLine(gl, stsColor, highlighted, rotatedPoints, zIndex);
    }

    public void mouseSelectedEdit(StsMouse mouse)
    {
        logMessage();
    }

    public void showPopupMenu(StsGLPanel glPanel, StsMouse mouse) { }

    public void pickVertices(StsGLPanel3d glPanel3d)
    {
        GL gl = glPanel3d.getGL();
        StsPoint[] rotatedPoints = computeRotatedCoorVertexPoints();
        for (int n = 0; n < rotatedPoints.length; n++)
        {
            gl.glInitNames();
            gl.glPushName(n);
            StsGLDraw.drawPoint(rotatedPoints[n].getXYZorT(), gl, 4);
            gl.glPopName();

        }
    }

    public void debugPrint()
    {
        System.out.println("StsLine mainDebug output. Number of vertices: " + lineVertices.getSize());

        try
        {
            StsSurfaceVertex vertex = (StsSurfaceVertex) lineVertices.getFirst();
            for (int nv = 0; vertex == null; nv++)
            {
                System.out.println("StsSurfaceVertex " + nv + " index " + vertex.getIndex());
                vertex = (StsSurfaceVertex) lineVertices.getNext();
            }
        }
        catch (Exception e)
        {
            System.out.println("StsSurfaceVertex.debugPrint() exception:\n" + e);
        }
    }

    /**
     * Methods for handling well highlighting.
     *
     * @param state indicates true or false
     * @return return true if toggled.
     */
    public boolean setHighlight(boolean state)
    {
        if (highlighted != state)
        {
            highlighted = state;
            return true;
        }
        else
        {
            return false;
        }
    }

    /** Draw all wells currently in highlighted list */
    static public void drawHighlightedLines(StsGLPanel3d glPanel3d)
    {
        //    	if(highlightedList != null) highlightedList.display(glPanel3d);
    }

    /** Clear the highlighted list */
    static public void clearHighlightedLines()
    {
        if (highlightedList != null)
        {
            highlightedList.clear();
        }
    }

    public boolean areAllSectionsVertical()
    {
        if (connectedSections == null)
        {
            return false;
        }
        for (int n = 0; n < connectedSections.getSize(); n++)
        {
            StsSection connectedSection = (StsSection) connectedSections.getElement(n);
            if (!connectedSection.isVertical())
            {
                return false;
            }
        }
        return true;
    }

    public boolean initialize(StsModel model)
    {
        //		return true;
        return initialize();
    }

    /**
     * Initialize well even if on uninitialized section. Return true only if
     * not on section or section is initialized
     */
    public boolean initialize()
    {
        if (initialized) return true;
        if (!computePoints()) return false;
        initialized = initializeSection();
        
        if(lineVertices != null)
            numberOfElements = lineVertices.getSize();

		StsSeismicVelocityModel velocityModel = currentModel.getProject().getSeismicVelocityModel();
        if (velocityModel != null) adjustFromVelocityModel(velocityModel);

        return initialized;
    }

    public boolean initializeSection()
    {
        if (onSection == null) return true;
        if (!onSection.initialized) return false;
        if (sectionInitialized) return true;
        projectToSection();
        sectionInitialized = true;
        return true;
    }

    public boolean reinitialize()
    {
        try
        {
            if (onSection == null)
            {
                return true;
            }
            projectToSection();
            if (onSection.initialized && !sectionInitialized)
            {
                sectionInitialized = true;
            }
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsLine.reinitialize() failed.",
                e, StsException.WARNING);
            return false;
        }
    }

    public void removeSection(StsSection section)
    {
        if (connectedSections != null)
        {
            deleteConnectedSection(section);
        }
    }

    static public StsLine buildVertical(StsSurfaceVertex vertex, byte type)
    {
        StsLine line;

        line = vertex.getSectionLine();
        if (line != null)
        {
            return line;
        }
        StsGridPoint gridPoint = new StsGridPoint(vertex.getPoint(), null);
        line = buildVertical(gridPoint, type);
        vertex.setAssociation(line);
        return line;
    }

    static public StsLine buildVertical(StsPoint point, byte type)
    {
        StsGridPoint gridPoint = new StsGridPoint(point, null);
        StsLine line = buildVertical(gridPoint, type);
        return line;
    }

    /*
        static public StsLine buildVertical(StsPoint point, boolean persistent)
        {
            StsLine line = new StsLine(persistent);
            StsGridPoint gridPoint = new StsGridPoint(point, null);
            line.constructVertical(gridPoint);
            return line;
        }
    */
    static public StsLine buildVertical(StsGridPoint gridPoint, byte type)
    {
        try
        {
            StsLine line = new StsLine(false);
            line.setZDomainOriginal(currentModel.getProject().getZDomain());
            line.constructVertical(gridPoint, type);
            line.addToModel();
            return line;
        }
        catch (Exception e)
        {
            StsException.systemError("StsLine.buildVertical(gridPoint) failed.");
            return null;
        }
    }

    /*
       static public StsLine buildVertical(StsSurfaceVertex vertex, byte type)
       {
           StsLine line;

           line = vertex.getSectionLine();
           if(line != null) return line;
     StsGridPoint gridPoint = new StsGridPoint(vertex.getPoint(), null);
           line = buildVertical(gridPoint);
           vertex.setAssociation(line);
           return line;
       }

       static public StsLine buildVertical(StsPoint point)
       {
          return buildVertical(point, true);
       }

       static public StsLine buildVertical(StsPoint point, boolean persistent)
       {
           StsGridPoint gridPoint = new StsGridPoint(point, null);
           StsLine line = buildVertical(gridPoint, persistent);
           return line;
       }

       static public StsLine buildVertical(StsGridPoint gridPoint)
       {
           return buildVertical(gridPoint, true);
       }

       static public StsLine buildVertical(StsGridPoint gridPoint, boolean persistent)
       {
           try
           {
               StsLine line = new StsLine(persistent);
               line.constructVertical(gridPoint, persistent);
               return line;
           }
           catch(Exception e)
           {
               StsException.systemError("StsLine.buildVertical(gridPoint) failed.");
               return null;
           }
       }
    */
    /** construct a vertical sectionLine thru this vertex */

    protected boolean constructVertical(StsGridPoint gridPoint, byte type)
    {
        return constructVertical(gridPoint, type, true);
    }

    /*
        private boolean constructVertical(StsGridPoint gridPoint, byte type, boolean persistent)
        {
            float z;
            if (lineVertices == null)
                lineVertices = StsObjectRefList.constructor(2, 1, "lineVertices", this, persistent);

            points = new StsPoint[2];
            StsProject project = currentModel.getProject();

            float x = gridPoint.getPoint().getX();
            float y = gridPoint.getPoint().getY();

            z = project.getZorTMin();
            points[0] = new StsPoint(x, y, z);
            addLineVertex(points[0], true, false);
            StsGridSectionPoint gridSectionPoint = new StsGridSectionPoint(points[0], row, col, null, null, persistent);
            vertex = new StsSurfaceVertex(gridSectionPoint, this, persistent); // point is cloned
            z = project.getZorTMax();
            points[1] = new StsPoint(x, y, z);
            return true;
        }
    */
    private boolean constructVertical(StsGridPoint gridPoint, byte type, boolean persistent)
    {
        StsPoint point;
        StsGridSectionPoint gridSectionPoint;
        StsSurfaceVertex vertex;

        this.type = type;

        if (lineVertices == null)
        {
            lineVertices = StsObjectRefList.constructor(2, 1, "lineVertices", this, persistent);
        }
        rotatedPoints = new StsPoint[2];
        int row = gridPoint.row;
        int col = gridPoint.col;

        StsProject project = currentModel.getProject();
        float[] rotatedXY = gridPoint.getPoint().v;

        zDomainSupported = project.getZDomainSupported();

        float t = 0.0f;
        float d = 0.0f;
        StsSeismicVelocityModel velocityModel = project.getSeismicVelocityModel();
        if (velocityModel != null)
        {
            t = velocityModel.getTimeMin();
            d = (float) velocityModel.getZ(rotatedXY[0], rotatedXY[1], t);
        }
        else
        {
            t = project.getTimeMin();
            d = project.getDepthMin();
        }

        //        points[0] = createStsPoint(rotatedXYZ[0], rotatedXYZ[1], t, d);
        //        point = createStsPoint(unrotatedXY[0], unrotatedXY[1], t, d);
        point = createStsPoint(rotatedXY[0], rotatedXY[1], t, d);
        rotatedPoints[0] = point;
        gridSectionPoint = new StsGridSectionPoint(point, row, col, null, null, persistent);
        vertex = new StsSurfaceVertex(gridSectionPoint, this, persistent); // point is cloned
        lineVertices.add(vertex);

        if (velocityModel != null)
        {
            t = velocityModel.getTimeMax();
            d = (float) velocityModel.getZ(rotatedXY[0], rotatedXY[1], t);
        }
        else
        {
            t = project.getTimeMax();
            d = project.getDepthMax();
        }
        //        points[1] = createStsPoint(rotatedXYZ[0], rotatedXYZ[1], t, d);
        //        point = createStsPoint(unrotatedXY[0], unrotatedXY[1], t, d);
        point = createStsPoint(rotatedXY[0], rotatedXY[1], t, d);
        rotatedPoints[1] = point;
        gridSectionPoint = new StsGridSectionPoint(point, row, col, null, null, persistent);
        vertex = new StsSurfaceVertex(gridSectionPoint, this, persistent); // point is cloned
        lineVertices.add(vertex);

        isVertical = true;
        if (stsColor == null)
        {
            stsColor = currentModel.getSpectrumClass().getCurrentSpectrumColor("Basic");
        }
        numberOfElements = lineVertices.getSize();

        return true;
    }

    /*
        public StsPoint createStsPoint(float x, float y, float z)
        {
            if (isDepth)
            {
                return new StsPoint(x, y, z);
            }
            else
            {
                StsPoint point = new StsPoint(5);
                point.setX(x);
                point.setY(y);
                point.setZ(z); // hack: we need to use a velocity, td curve, or td model
                point.setT(z);
                return point;
            }
        }
    */
    public StsPoint createStsPoint(float x, float y, float t, float d)
    {
        /*
                if (isDepth)
                {
                    return new StsPoint(x, y, z);
                }
                else
                {
        */
        StsPoint point = new StsPoint(5);
        point.setX(x);
        point.setY(y);
        point.setZ(d);
        point.setT(t);
        return point;
        //         }
    }

    public StsPoint createStsPoint(StsPoint point)
    {
        return createStsPoint(point.v);
    }

    public StsPoint createStsPoint(float[] xyz)
    {
        if (isDepth)
        {
            return new StsPoint(xyz[0], xyz[1], xyz[2]);
        }
        else
        {
            StsPoint point = new StsPoint(5);
            point.setX(xyz[0]);
            point.setY(xyz[1]);
            point.setZ(xyz[2]);
            if (xyz.length < 5)
            {
                point.setT(xyz[2]);
            }
            else
            {
                point.setT(xyz[4]);
            }
            return point;
        }
    }

    public boolean construct(StsSurfaceVertex[] vertices)
    {
        if (vertices == null) return false;
        int nVertices = vertices.length;
        if (lineVertices == null)
        {
            lineVertices = StsObjectRefList.constructor(nVertices, 2, "lineVertices", this);

        }
        lineVertices.add(vertices);
        numberOfElements = lineVertices.getSize();
        zDomainSupported = currentModel.getProject().getZDomainSupported();
        return true;
    }

    /*
        public boolean construct(StsSurfaceVertex[] vertices)
        {
            StsPoint point;
            StsGridSectionPoint gridSectionPoint;
            StsSurfaceVertex vertex;

            if (vertices == null) return false;

            int nVertices = vertices.length;
            StsGridSectionPoint[] gridPoints = new StsGridSectionPoint[nVertices];
            for(int n = 0; n < nVertices; n++)
                    gridPoints[n] = vertices[n].getSurfacePoint();

            int nPoints = gridPoints.length;
            if (nPoints == 0)
            {
                return false;
            }
            if (lineVertices == null)
            {
                lineVertices = StsObjectRefList.constructor(nPoints, 1, "lineVertices", this);
            }
            points = new StsPoint[nPoints];
            StsProject project = currentModel.getProject();
            xOrigin = project.getXOrigin();
            yOrigin = project.getYOrigin();
            for (int n = 0; n < nPoints; n++)
            {
                float rowF = gridPoints[n].getGridRowF();
                float colF = gridPoints[n].getGridColF();

                point = gridPoints[n].point.copy(); // make a copy of rotated point; will subsequently unrotate it
    //            gridPoints[n].point = point;
                points[n] = point;
                float[] unrotatedXY = getUnrotatedRelativeXYFromRotatedXY(point.getX(), point.getY());
    //            point = gridPoints[n].point.copy();
                point.setX(unrotatedXY[0]);
                point.setY(unrotatedXY[1]);
                gridSectionPoint = new StsGridSectionPoint(point, rowF, colF, null, null, true);
                vertex = new StsSurfaceVertex(gridSectionPoint, this); // point is cloned
                lineVertices.add(vertex);
                if (stsColor == null)
                {
                    stsColor = currentModel.getSpectrumClass().getCurrentSpectrumColor("Basic");
                }
            }
            return true;
        }
    */
    private float[] getUnrotatedRelativeXY(float x, float y)
    {
        return currentModel.getProject().getUnrotatedRelativeXYFromRotatedXY(x, y);
    }

    public StsSurfaceVertex getConstructSurfaceEdgeVertex(StsSurface surface, StsBlock block)
    {
        StsSurfaceVertex vertex;

        vertex = getSurfaceBlockVertex(surface, block);
        if (vertex != null) return vertex;

        // if not, copy point from an existing sectionLineVertex
        StsPoint sectionVertexPoint;
        StsSurfaceVertex sectionVertex = getSurfaceEdgeVertex(surface);
        if (sectionVertex != null)
            sectionVertexPoint = sectionVertex.getPoint();
        else
        {
            StsGridPoint gridPoint = computeGridIntersect(surface);
            if (gridPoint == null)
            {
                StsException.systemError("StsLine.construtInitialEdgeVertex() failed for " + getName() + " intersecting surface " + surface.getName());
                gridPoint = computeGridIntersect(surface);
                return null;
            }
            sectionVertexPoint = gridPoint.getPoint();
        }
//         point = sectionVertexPoint.getXYZorTPoint();
        vertex = new StsSurfaceVertex(sectionVertexPoint, this, surface, block, false);
        addSurfaceEdgeVertex(vertex);
        return vertex;
    }

    public void addSectionEdgeVertex(StsSurfaceVertex vertex)
    {
        if(sectionEdgeVertices == null) sectionEdgeVertices = StsObjectRefList.constructor(2, 2, "sectionEdgeVertices", this);
        sectionEdgeVertices.add(vertex);
    }

    public void addSurfaceEdgeVertex(StsSurfaceVertex vertex)
    {
        if(surfaceEdgeVertices == null) surfaceEdgeVertices = new StsObjectList(2, 2);
        surfaceEdgeVertices.add(vertex);

    }
    /*
        public StsSurfaceVertex constructInitialEdgeVertex(StsModelSurface surface, StsBlock block)
        {
            StsSurfaceVertex vertex;
            StsPoint point;

            vertex = getSurfaceBlockVertex(surface, block);
            if(vertex != null) return vertex;

            // if not, copy from an existing sectionLineVertex
            point = getSurfaceVertex(surface).getPoint();
            vertex = new StsSurfaceVertex(point, this, surface, block);
            surfaceVertices.add(vertex);
            return vertex;
        }
    */
    private StsSurfaceVertex getSurfaceBlockVertex(StsSurface surface, StsBlock block)
    {
        if (surfaceEdgeVertices == null) return null;
        int nVertices = surfaceEdgeVertices.getSize();
        for (int n = 0; n < nVertices; n++)
        {
            StsSurfaceVertex vertex = (StsSurfaceVertex) surfaceEdgeVertices.getElement(n);
            if (vertex.getSurface() == surface && vertex.getBlock() == block) return vertex;
        }
        return null;
    }

    public StsSurfaceVertex getVertexBelow(StsSurfaceVertex vertex)
    {
        StsXYSurfaceGridable surface = vertex.getSurface();
        if (!(surface instanceof StsModelSurface))
        {
            return null;
        }
        StsModelSurface modelSurface = (StsModelSurface) surface;
        StsSurface surfaceBelow = modelSurface.getModelSurfaceBelow();
        if (surfaceBelow == null)
        {
            return null;
        }

        StsBlock block = vertex.getBlock();

        int nVertices = lineVertices.getSize();
        for (int n = 0; n < nVertices; n++)
        {
            vertex = (StsSurfaceVertex) lineVertices.getElement(n);
            if (vertex.getSurface() == surfaceBelow && vertex.getBlock() == block)
            {
                return vertex;
            }
        }
        return null;
    }

    public StsLineZone getLineZone(StsSurfaceVertex topVertex, StsSurfaceVertex botVertex,
                                   StsSection onSection, int sectionSide, int direction)
    {
        if (onSection == null)
        {
            return null;
        }

        // If this is a defining well for the section, return wellZone on sectionSide
        if (onSection.getFirstLine() == this || onSection.getLastLine() == this)
        {
            return getLineZone(topVertex, botVertex, onSection, sectionSide);
        }

        int position = getLineSectionEnd();
        if (position == NONE)
        {
            return null;
        }

        if (position == MINUS && direction == MINUS || position == PLUS && direction == PLUS)
        {
            return getLineZone(topVertex, botVertex, onSection, RIGHT);
        }
        else
        {
            return getLineZone(topVertex, botVertex, onSection, RIGHT);
        }
    }

    public StsLineZone getLineZone(StsSurfaceVertex topVertex, StsSurfaceVertex botVertex, StsSection onSection,
                                   int sectionSide)
    {
        StsLineZone zone;

        if (zones == null)
        {
            zones = new StsList(2, 2);
        }

        int nZones = zones.getSize();
        for (int n = 0; n < nZones; n++)
        {
            zone = (StsLineZone) zones.getElement(n);
            if (zone.getTop() == topVertex)
            {
                return zone;
            }
        }

        /** Couldn't find zone: make one and store it in zones */
        zone = new StsLineZone(currentModel, topVertex, botVertex, this);
        zones.add(zone);
        return zone;
    }

    public void deleteTransients()
    {
        zones = null;
    }

    public float getDipAngle(float[] dipDirectionVector, StsPoint linePoint, int sectionEnd)
    {
        StsPoint topPoint, botPoint, lineVector;

        if (onSection != null)
        {
            StsPoint tangent = onSection.getTangentAtPoint(linePoint, sectionSide);
            if (sectionEnd == MINUS && sectionSide == RIGHT || sectionEnd == PLUS && sectionSide == LEFT)
            {
                tangent.reverse();
            }
            dipDirectionVector = tangent.v;
        }

        topPoint = getTopPoint();
        botPoint = getBotPoint();

        lineVector = new StsPoint(botPoint);
        lineVector.subtract(topPoint);
        lineVector.normalize();

        float horizontal = StsMath.dot(dipDirectionVector, lineVector.v);
        float vertical = lineVector.getZorT();
        float dipAngle = StsMath.atan2(vertical, horizontal);
        if (dipAngle > 180.0f)
        {
            dipAngle -= 360.0f; // reset dipAngle between -180 to +180
        }
        return dipAngle;
    }

    public void adjustDipAngle(float dipAngleChange, float[] axis, StsPoint linePoint)
    {
        StsRotationMatrix rotMatrix;

        rotMatrix = StsRotationMatrix.constructRotationMatrix(linePoint.v, axis, dipAngleChange);

        int nVertices = lineVertices.getSize();
        for (int n = 0; n < nVertices; n++)
        {
            StsSurfaceVertex vertex = (StsSurfaceVertex) lineVertices.getElement(n);
            StsPoint point = vertex.getPoint();
            rotMatrix.pointRotate(point.v);
        }
        computePointsProjectToSection();
        setIsVertical(false);
    }

    public StsPoint getPointOnLineNearestMouse(int nSegment, StsMouse mouse, StsGLPanel3d glPanel3d)
    {
        if (rotatedPoints == null)
        {
            return null;
        }
        int nPoints = rotatedPoints.length;
        int n = Math.min(nPoints - 1, nSegment);
        return glPanel3d.getPointOnLineNearestMouse(mouse, rotatedPoints[n], rotatedPoints[n + 1]);
    }

    public boolean extendEnds()
    {
        float z;
        StsPoint point;
        float zTop = this.getTopLineVertex().getPoint().getZorT();
        float zBot = this.getBotLineVertex().getPoint().getZorT();
        float zDif = zBot - zTop;
        zTop = zTop - zDif / 2;
        zBot = zBot + zDif / 2;
        StsProject project = currentModel.getProject();

        int sectionIndexMin = project.getIndexAbove(zTop);
        int sectionIndexMax = project.getIndexBelow(zBot);
        zTop = project.getZAtIndex(sectionIndexMin);
        zBot = project.getZAtIndex(sectionIndexMax);

        return computeExtendedPoints(zTop, zBot);
    }

    /*
// Extend ends to top and bottom of project
    public void extendEnds()
    {
     float z;
     StsPoint point;

           computePoints();
     StsProject project = currentModel.getProject();
     z = project.getZMin();
     point = getPointAtZ(z, true);
           addLineVertex(point, true);
     z = project.getZMax();
     point = getPointAtZ(z, true);
           addLineVertex(point, true);
           computePoints();
    }
    */
    public StsBoundingBox getBoundingBox()
    {
        if (boundingBox != null)
        {
            return boundingBox;
        }
        return new StsBoundingBox(rotatedPoints, xOrigin, yOrigin);
    }

    public void addToRotatedBoundingBox(StsRotatedBoundingBox rotatedBoundingBox)
    {
        rotatedBoundingBox.addPoints(rotatedPoints, xOrigin, yOrigin);
    }

    public void logMessage()
    {
        logMessage(lineOnSectionLabel());
    }

    /**
     * set x-y-z geometry and z-measured depth pairs
     *
     * @param logVectors set of logVectors including x, y, depth, and mDepth
     * @return true if executes correctly
     *         <p/>
     *         If project origin hasn't been set, use origin of this well as project origin;
     *         x and y vectors need not be adjusted.  If origin has been set, then compute
     *         offset of this well's origin from project origin and add to x & y vectors.
     */
    public boolean setVectors(StsLogVector[] logVectors)
    {
        StsLogVector xLogVector, yLogVector, zLogVector, mLogVector;
        float[] xArray, yArray, zArray, mArray;

        checkExcursion = false;
        //useDataStore = true;
        try
        {
            // get project Z tolerance to limit number of lineVertices we create
            //            StsProject project = currentModel.getProject();

            // set X-Y-Z points
            xLogVector = StsLogVector.getVectorOfType(logVectors, StsLogVector.X);
            if (xLogVector == null)
            {
                return false;
            }
            xArray = xLogVector.getFloats();
            if (xArray == null)
            {
                return false;
            }

            yLogVector = StsLogVector.getVectorOfType(logVectors, StsLogVector.Y);
            if (yLogVector == null)
            {
                return false;
            }
            yArray = yLogVector.getFloats();
            if (yArray == null)
            {
                return false;
            }

            zLogVector = StsLogVector.getVectorOfType(logVectors, StsLogVector.DEPTH);
            if (zLogVector == null)
            {
                return false;
            }
            zArray = zLogVector.getFloats();
            if (zArray == null)
            {
                return false;
            }

            mLogVector = StsLogVector.getVectorOfType(logVectors, StsLogVector.MDEPTH);
            if (mLogVector == null)
            {
                return false;
            }
            mArray = mLogVector.getFloats();
            if (zArray == null)
            {
                return false;
            }
            int nValues = xArray.length;

            xOrigin = xLogVector.getOrigin();
            yOrigin = yLogVector.getOrigin();
            /*
              StsBoundingBox unrotatedBoundingBox = project.getUnrotatedBoundingBox();
                       if(!unrotatedBoundingBox.originSet)
                           unrotatedBoundingBox.setOrigin(xOrigin, yOrigin);
                       else
                       {
                           xLogVector.checkAdjustOrigin(unrotatedBoundingBox.getXOrigin());
                           yLogVector.checkAdjustOrigin(unrotatedBoundingBox.getYOrigin());
                       }
            */
            // set bounding box x and y range
            /*
               project.adjustBoundingBoxXRange(xLogVector.getMinValue(), xLogVector.getMaxValue());
               project.adjustBoundingBoxYRange(yLogVector.getMinValue(), yLogVector.getMaxValue());
               project.adjustBoundingBoxZRange(zLogVector.getMinValue(), zLogVector.getMaxValue());
            */
            //            project.adjustRange(getRelativeRange(xLogVector, yLogVector, zLogVector));

            //            zInc = project.getZInc();
            if (lineVertices == null)
            {
                lineVertices = StsObjectRefList.constructor(nValues, 1, "lineVertices", this);
            }
            StsPoint point = new StsPoint(xArray[0], yArray[0], zArray[0], mArray[0]);
            lineVertices.add(new StsSurfaceVertex(point, this));
            int nvm1 = nValues - 1;
            for (int i = 1; i < nvm1; i++)
            {
                point.setX(xArray[i]);
                point.setY(yArray[i]);
                point.setZ(zArray[i]);
                point.setF(mArray[i]);
                lineVertices.add(new StsSurfaceVertex(point, this));
            }
            point.setX(xArray[nvm1]);
            point.setY(yArray[nvm1]);
            point.setZ(zArray[nvm1]);
            point.setF(mArray[nvm1]);
            lineVertices.add(new StsSurfaceVertex(point, this));

            /** If we don't have a color yet, get it now */
            if (stsColor == null)
            {
                stsColor = currentModel.getCurrentSpectrumColor("Basic");
            }
        }
        catch (Exception e)
        {
            StsMessageFiles.logMessage("Unable to build well line geometry");
            return false;
        }
        return true;
    }

    public float[] getRelativeRange(StsLogVector xLogVector, StsLogVector yLogVector, StsLogVector zLogVector)
    {
        float xMin = xLogVector.getMinValue();
        float xMax = xLogVector.getMaxValue();
        float yMin = yLogVector.getMinValue();
        float yMax = yLogVector.getMaxValue();
        float zMin = zLogVector.getMinValue();
        float zMax = zLogVector.getMaxValue();

        return new float[]
            {
                xMin, xMax, yMin, yMax, zMin, zMax};
    }

    public boolean addToProject(byte zDomainSupported)
    {
        this.zDomainSupported = zDomainSupported;
        return currentModel.getProject().addToProject(this);
        //		project.checkSetOrigin(xOrigin, yOrigin);
        //		StsBoundingBox unrotatedBoundingBox = project.getUnrotatedBoundingBox();
        //        addToUnrotatedBoundingBox(unrotatedBoundingBox);
        /*
          int nVertices = lineVertices.getSize();
          boolean isDepth = project.isDepth();
          for(int n = 0; n < nVertices; n++)
          {
           StsSurfaceVertex vertex = (StsSurfaceVertex)lineVertices.getElement(n);
           StsPoint point = vertex.getPoint();
           if(isDepth)
            unrotatedBoundingBox.addPoint(point, xOrigin, yOrigin);
           else
            unrotatedBoundingBox.addPoint(point.getX(), point.getY(), point.getT(), xOrigin, yOrigin);
          }
        */
    }

    public StsFieldBean[] getDisplayFields()
    {
        if (pseudoDisplayFields == null)
        {
            pseudoDisplayFields = new StsFieldBean[]
                {
                    new StsBooleanFieldBean(StsLine.class, "isVisible", "Enable"),
                    new StsBooleanFieldBean(StsLine.class, "drawZones", "Zones"),
                    new StsColorComboBoxFieldBean(StsLine.class, "stsColor", "Color"),
                    new StsFloatFieldBean(StsLine.class, "topZ", false, "Min Depth"),
                    new StsFloatFieldBean(StsLine.class, "botZ", false, "Max Depth"),
                    new StsDoubleFieldBean(StsLine.class, "xOrigin", false, "X Origin"),
                    new StsDoubleFieldBean(StsLine.class, "yOrigin", false, "Y Origin")
                };
        }
        return pseudoDisplayFields;
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

    static public StsFieldBean[] getStaticDisplayFields()
    {
        return pseudoDisplayFields;
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
        currentModel.getCreateStsClass("com.Sts.DBTypes.StsLine").selected(this);
        this.setStsColor(StsColor.BLACK);
    }

    public StsPoint[] getTimePoints(StsPoint[] depthPoints)
    {
        return null;
    }

    public float[] getDepthArrayFromVertices()
    {
        if (lineVertices == null)
        {
            return null;
        }
        int nVertices = lineVertices.getSize();
        float[] depths = new float[nVertices];
        for (int n = 0; n < nVertices; n++)
        {
            StsSurfaceVertex vertex = (StsSurfaceVertex) lineVertices.getElement(n);
            StsPoint point = vertex.getPoint();
            depths[n] = point.getZ();
        }
        return depths;
    }

    public float[] getMDepthArrayFromVertices()
    {
        if (lineVertices == null)
        {
            return null;
        }
        int nVertices = lineVertices.getSize();
        float[] mdepths = new float[nVertices];
        for (int n = 0; n < nVertices; n++)
        {
            StsSurfaceVertex vertex = (StsSurfaceVertex) lineVertices.getElement(n);
            StsPoint point = vertex.getPoint();
            mdepths[n] = point.getM();
        }
        return mdepths;
    }

    static public void addMDepthToPoints(StsPoint[] points)
    {
        if (points == null) return;
        StsPoint point1 = points[0];
        float mdepth = 0.0f;
        point1.setM(0.0f);
        for (int n = 1; n < points.length; n++)
        {
            StsPoint point0 = point1;
            point1 = points[n];
            mdepth += point0.distance(point1);
            point1.setM(mdepth);
        }
    }

    static public void addZorTasMDepthToPoints(StsPoint[] points)
    {
        if (points == null) return;
        for (int n = 0; n < points.length; n++)
            points[n].setM(points[n].getZorT());
    }

    static public void addMDepthToVertices(StsSurfaceVertex[] vertices)
    {
        if (vertices == null) return;
        StsPoint point1 = vertices[0].getPoint();
        float mdepth = 0.0f;
        point1.setM(0.0f);
        for (int n = 1; n < vertices.length; n++)
        {
            StsPoint point0 = point1;
            point1 = vertices[n].getPoint();
            mdepth += point0.distance(point1);
            point1.setM(mdepth);
        }
    }

    public ArrayList getGridCrossingPoints(StsXYSurfaceGridable grid)
    {
        int nPoints = rotatedPoints.length;
        StsGridCrossingPoint[] gridPoints = new StsGridCrossingPoint[nPoints];
        for (int n = 0; n < nPoints; n++)
        {
            gridPoints[n] = new StsGridCrossingPoint(grid, rotatedPoints[n]);
        }
        return getGridCrossingPoints(gridPoints);
    }

    static public ArrayList getGridCrossingPoints(StsGridCrossingPoint[] gridPoints)
    {
        ArrayList list = new ArrayList();
        int nPoints = gridPoints.length;
        StsGridCrossingPoint gridPoint1 = gridPoints[0];
        list.add(gridPoint1);
        StsGridCrossingPoint lastGridPoint = gridPoint1;
        for (int n = 1; n < nPoints; n++)
        {
            StsGridCrossingPoint gridPoint0 = gridPoint1;
            gridPoint1 = gridPoints[n];
            StsGridCrossings gridCrossings = new StsGridCrossings(gridPoint0, gridPoint1, true);
            ArrayList gridCrossingPoints = gridCrossings.gridPoints;
            int nCrossings = gridCrossingPoints.size();
            for (int i = 0; i < nCrossings; i++)
            {
                StsGridCrossingPoint gridPoint = (StsGridCrossingPoint) gridCrossingPoints.get(i);
                checkAddGridPoint(gridPoint, lastGridPoint, list);
                lastGridPoint = gridPoint;
            }
        }
        if (gridPoint1.rowOrCol == StsParameters.NONE)
        {
            list.add(gridPoint1);
        }
        return list;
    }

    static private void checkAddGridPoint(StsGridCrossingPoint gridPoint, StsGridCrossingPoint lastGridPoint, ArrayList list)
    {
        if (gridPoint.rowOrCol == StsParameters.NONE)
        {
            return;
        }
        if (gridPoint.sameAs(lastGridPoint))
        {
            return;
        }
        list.add(gridPoint);
    }

    static public ArrayList getCellCrossingPoints(StsXYGridable grid, StsPoint[] points)
    {
        //        StsPoint[] points = getRotatedPoints();
        int nPoints = points.length;
        StsGridCrossingPoint[] gridPoints = new StsGridCrossingPoint[nPoints + 2];
        // extend the curtain laterally at the top and bottom so the viewer has good coverage
        // if the well is vertical, just make it 20 traces wide in inline direction
        for (int n = 1; n <= nPoints; n++)
        {
            gridPoints[n] = new StsGridCrossingPoint(grid, points[n - 1]);
            gridPoints[n].adjustToCellGridding();
        }
        int maxNRow = grid.getNRows() - 1;
        int maxNCol = grid.getNCols() - 1;
        float dRowF, dColF;
        float rowF1, rowF2, colF1, colF2;
        rowF1 = gridPoints[1].iF;
        rowF2 = gridPoints[nPoints].iF;
        dRowF = rowF2 - rowF1;
        colF1 = gridPoints[1].jF;
        colF2 = gridPoints[nPoints].jF;
        dColF = colF2 - colF1;
        // well is vertical: extend in inline direction
        if (Math.abs(dRowF) < 0.01f && Math.abs(dColF) < 0.01f)
        {
            rowF1 = Math.max(0, rowF1 - 10);
            gridPoints[0] = new StsGridCrossingPoint(rowF1, colF1, grid);
            rowF2 = Math.min(maxNRow, rowF1 + 20);
            gridPoints[nPoints + 1] = new StsGridCrossingPoint(rowF2, colF2, grid);
        }
        // well is almost vertical: extend in direction of plane thru top and bottom points
        else if (Math.abs(dRowF) < 1.0f && Math.abs(dColF) < 1.0f || nPoints <= 4)
        {
            float d = (float) Math.sqrt(dRowF * dRowF + dColF * dColF);
            float ratio = 10 / d;

            rowF1 -= dRowF * ratio;
            rowF1 = StsMath.minMax(rowF1, 0, maxNRow);
            colF1 -= dColF * ratio;
            colF1 = StsMath.minMax(colF1, 0, maxNCol);
            gridPoints[0] = new StsGridCrossingPoint(rowF1, colF1, grid);

            rowF2 += dRowF * ratio;
            rowF2 = StsMath.minMax(rowF2, 0, maxNRow);
            colF2 += dColF * ratio;
            colF2 = StsMath.minMax(colF2, 0, maxNCol);
            gridPoints[nPoints + 1] = new StsGridCrossingPoint(rowF2, colF2, grid);
        }
        else // search up from bottom and top until we have enough offset to locate an added point
        {
            rowF1 = gridPoints[1].iF;
            colF1 = gridPoints[1].jF;
            for (int n = 2; n < nPoints / 2; n++)
            {
                rowF2 = gridPoints[n].iF;
                dRowF = rowF2 - rowF1;
                colF2 = gridPoints[n].jF;
                dColF = colF2 - colF1;
                if (Math.abs(dRowF) > 1.0f || Math.abs(dColF) > 1.0f)
                {
                    break;
                }
            }
            float d = (float) Math.sqrt(dRowF * dRowF + dColF * dColF);
            float ratio = 10 / d;

            rowF1 -= dRowF * ratio;
            rowF1 = StsMath.minMax(rowF1, 0, maxNRow);
            colF1 -= dColF * ratio;
            colF1 = StsMath.minMax(colF1, 0, maxNCol);
            gridPoints[0] = new StsGridCrossingPoint(rowF1, colF1, grid);

            rowF1 = gridPoints[nPoints].iF;
            colF1 = gridPoints[nPoints].jF;
            for (int n = nPoints - 1; n >= nPoints / 2; n--)
            {
                rowF2 = gridPoints[n].iF;
                dRowF = rowF2 - rowF1;
                colF2 = gridPoints[n].jF;
                dColF = colF2 - colF1;
                if (Math.abs(dRowF) > 1.0f || Math.abs(dColF) > 1.0f)
                {
                    break;
                }
            }
            d = (float) Math.sqrt(dRowF * dRowF + dColF * dColF);
            ratio = 10 / d;

            rowF1 -= dRowF * ratio;
            rowF1 = StsMath.minMax(rowF1, 0, maxNRow);
            colF1 -= dColF * ratio;
            colF1 = StsMath.minMax(colF1, 0, maxNCol);
            gridPoints[nPoints + 1] = new StsGridCrossingPoint(rowF1, colF1, grid);
        }
        return getGridCrossingPoints(gridPoints);
    }

    /*
         public boolean adjustTimePoints(StsSeismicVelocityModel velocityModel)
         {
             StsPoint[] points = getRotatedPoints();

             if (points == null)
             {
                 computeXYZPoints();
                 points = getRotatedPoints();
             }
             if (points == null) return false;

             for (int n = 0; n < points.length; n++)
             {
                 try
                 {
                     float t = (float)velocityModel.getT(points[n].v);
                     points[n].setT(t);
                 }
                 catch (Exception e)
                 {
                     StsMessageFiles.errorMessage("Failed to adjust well points for well " + getName() + " points probably not in time.");
                     return false;
                 }

             }
             return true;
         }
    */
    private boolean adjustTimeOrDepthPoints(StsSeismicVelocityModel velocityModel, boolean isOriginalDepth)
    {
        StsPoint[] points = getRotatedPoints();

        if (points == null)
        {
            computePoints();
            points = getRotatedPoints();
        }
        if (points == null)
        {
            return false;
        }
        boolean success = true;
        for (int n = 0; n < points.length; n++)
            if (!velocityModel.adjustTimeOrDepthPoint(points[n], isOriginalDepth)) success = false;
        return success;
    }

    public boolean adjustTimeOrDepth(StsSeismicVelocityModel velocityModel)
    {
        return adjustTimeOrDepth(velocityModel, isOriginalDepth());
    }

    public boolean isOriginalDepth()
    {
        return getZDomainOriginal() == StsProject.TD_DEPTH;
    }

    public boolean adjustTimeOrDepth(StsSeismicVelocityModel velocityModel, boolean isOriginalDepth)
    {
        int nVertices = lineVertices.getSize();
        boolean success = true;
        for (int n = 0; n < nVertices; n++)
        {
            StsSurfaceVertex vertex = (StsSurfaceVertex) lineVertices.getElement(n);
            StsPoint point = vertex.getPoint();
            if (!velocityModel.adjustTimeOrDepthPoint(point, isOriginalDepth)) success = false;
        }
        float zMin = currentModel.getProject().getZorTMin();
        float zMax = currentModel.getProject().getZorTMax();
        this.computeExtendedPointsProjectToSection(zMin, zMax);
        adjustTimeOrDepthPoints(velocityModel, isOriginalDepth);
        return success;
    }

    public boolean isVerticesRotated()
    {
        return isVerticesRotated;
    }

    public void setVerticesRotated(boolean verticesRotated)
    {
        fieldChanged("isVerticesRotated", verticesRotated);
    }

    public byte getZDomainOriginal()
    {
        return zDomainOriginal;
    }

    public void setZDomainOriginal(byte zDomainOriginal)
    {
        this.zDomainOriginal = zDomainOriginal;
        this.zDomainSupported = zDomainOriginal;
    }

    public boolean checkAdjustFromVelocityModel()
    {
        StsSeismicVelocityModel velocityModel = currentModel.getProject().getSeismicVelocityModel();
        if (velocityModel == null) return false;
        return adjustFromVelocityModel(velocityModel);
    }
	/** adjust the lineVertices time or depth values using a new velocity model.
	 *  Message max time adjustment.
	 */
    public boolean adjustFromVelocityModel(StsSeismicVelocityModel velocityModel)
    {
        if (lineVertices == null) return false;
        try
        {
            StsProject project = currentModel.getProject();
            int nLineVertices = lineVertices.getSize();
            rotatedPoints = new StsPoint[nLineVertices];

            float dXOrigin = (float) (xOrigin - project.getXOrigin());
            float dYOrigin = (float) (yOrigin - project.getYOrigin());

            float maxChange = 0.0f;
            float[] times = new float[nLineVertices];
            for (int n = 0; n < nLineVertices; n++)
            {
                StsSurfaceVertex lineVertex = (StsSurfaceVertex) lineVertices.getElement(n);
                StsPoint originalVertexPoint = lineVertex.getPoint();
                StsPoint vertexPoint = new StsPoint(5, originalVertexPoint);
                if (!isVerticesRotated)
                {
                    float[] xy = project.getRotatedRelativeXYFromUnrotatedRelativeXY(dXOrigin + vertexPoint.v[0], dYOrigin + vertexPoint.v[1]);
                    vertexPoint.setX(xy[0]);
                    vertexPoint.setY(xy[1]);
                }
                if (zDomainOriginal == StsProject.TD_DEPTH)
                {
                    float t = (float) velocityModel.getT(vertexPoint.v);
                    float currentT = vertexPoint.getT();
                    maxChange = Math.max(maxChange, Math.abs(t - currentT));
                    vertexPoint.setT(t);
                    times[n] = t;
                    rotatedPoints[n] = vertexPoint;
					originalVertexPoint.setT(t);
                }
                else // zDomainOriginal is TIME
                {
                    float z = (float) velocityModel.getZ(vertexPoint);
                    vertexPoint.setZ(z);
                    rotatedPoints[n] = vertexPoint;
					originalVertexPoint.setZ(z);
                }
            }
            //             saveVertexTimesToDB(times);
            StsMessageFiles.logMessage("Well " + getName() + " max time adjustment " + maxChange);
            //             currentModel.addMethodCmd(this, "adjustTimePoints", new Object[] {velocityModel}, "adjustTimePoints for " + getName());
            return true;
        }
        catch (Exception e)
        {
            StsMessageFiles.errorMessage("Failed to adjust well " + getName() + " points probably not in time.");
            return false;
        }
    }

    public String toSectionString()
    {
        String sectionString = "none";
        if (onSection != null)
            sectionString = onSection.toString();
        return getName() + " section: " + sectionString;
    }

}
