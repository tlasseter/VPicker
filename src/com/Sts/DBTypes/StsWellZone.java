
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

public class StsWellZone extends StsObject implements StsSelectable, StsEdgeLinkable
{
    // instance fields
    protected StsWell well;
    protected StsSurfaceVertex topVertex = null;
    protected StsSurfaceVertex botVertex = null;
    protected StsObjectRefList edgePoints;
    protected float[] subZoneZs = null;
    protected boolean displayOnWellLine = true;
    protected StsWellZoneSet wellZoneSet;

    static final int NONE = StsParameters.NONE;

    /** constructors */

	public StsWellZone()
	{
	}

    public StsWellZone(StsWell well, int type, String name,
                       StsSurfaceVertex topVertex, StsSurfaceVertex botVertex)
    {
        if (name==null)
        {
            StsException.systemError("StsWellZone.StsWellZone:" +
                                     " Cannot create a well zone with a null name.");
        }
        if (well==null)
        {
            StsException.systemError("StsWellZone.StsWellZone:" +
                                     " Cannot create a zone for a null well.");
        }

	    this.well = well;
        wellZoneSet = getWellZoneSet(name, type);
        if (wellZoneSet==null) wellZoneSet = new StsWellZoneSet(type, name);
        wellZoneSet.addWellZone(this);
	    setTop(topVertex);
	    setBot(botVertex);
        well.addZone(this);
    }

    private StsWellZoneSet getWellZoneSet(String name, int type)
    {
        if (name == null) return null;
        StsClass zoneSetClass = currentModel.getCreateStsClass(StsWellZoneSet.class);
        if (zoneSetClass == null) return null;
        int nWellZoneSets = zoneSetClass.getSize();
        for (int i = 0; i < nWellZoneSets; i++)
        {
            StsWellZoneSet wellZoneSet = (StsWellZoneSet) zoneSetClass. getElement(i);
            if (name.equals(wellZoneSet.getName()) && type == wellZoneSet.getZoneType())
                return wellZoneSet;
        }
        return null;
    }

    /** constructor for 2 markers */
    public StsWellZone(StsWell well, int type, String name,
                       StsWellMarker topMarker, StsWellMarker baseMarker,
                       StsZone parentZone)
    {
        if (topMarker==null || baseMarker==null)
        {
            StsException.systemError("StsWellZone.StsWellZone:"
                      + " Top and/or base markers are null.");
        }
        if (well==null)
        {
            StsException.systemError("StsWellZone.StsWellZone:" +
                                     " Cannot create a zone for a null well.");
        }
	    this.well = well;

	    setTop(new StsSurfaceVertex(topMarker.getLocation(), null, null, null, true));
	    setBot(new StsSurfaceVertex(baseMarker.getLocation(), null, null, null, true));
//	    setTop(new StsSurfaceVertex(topMarker.getLocation(), well, null, null));
//	    setBot(new StsSurfaceVertex(baseMarker.getLocation(), well, null, null));
//        well.addZone(this);
    }

    static public StsWellZone constructor(StsWell well, int type,
                                          String name, StsWellMarker topMarker,
                                          StsWellMarker baseMarker, StsZone parentZone)
    {
        try{ return new StsWellZone(well, type, name, topMarker, baseMarker, parentZone); }
        catch(Exception e) {return null; }
    }

    // included for interface compatibility
    public int getRowOrCol() { return -1; }
    public int getRowCol() { return -1; }
    public StsEdgeLinkable getNextEdge() { return null; }
    public StsEdgeLinkable getPrevEdge() { return null; }
    public StsSurfaceVertex getNextVertex() { return botVertex; }
    public StsSurfaceVertex getPrevVertex() { return topVertex; }

    public boolean initialize(StsModel model) { return true; }

    // Accessors
    public StsWell getWell(){ return well; }
    public void setTop(StsSurfaceVertex topVertex)
    {
        if (topVertex == null)
        {
            StsException.systemError("StsWellZone.setTop:"
                      + " Cannot create a zone with a null vertex");
        }
        this.topVertex = topVertex;
        if (botVertex==null) return;

        float topZ = getTopZ();
        float baseZ = getBotZ();
        if (topZ >= baseZ)  // doesn't handle horizontal wells!
        {
            StsException.systemError("StsWellZone.setTop:"
                      + " Zone must have positive thickness.");
        }
        setSubZones();
//        well.setWellLineNeedsRebuild();
    }
    public StsSurfaceVertex getTop(){ return topVertex; }
    public float getTopZ()
    {
        if (topVertex==null) StsException.systemError("StsWellZone.getTopZ:"
                          + " Null top vertex.");
        return topVertex.getPoint().getZ();
    }
    public void setBot(StsSurfaceVertex botVertex)
    {
        if (botVertex == null)
        {
            StsException.systemError("StsWellZone.setBot:"
                      + " Cannot create a zone with a null vertex");
        }
        this.botVertex = botVertex;
        if (topVertex==null) return;

        float topZ = getTopZ();
        float baseZ = getBotZ();
        if (topZ >= baseZ)  // doesn't handle horizontal wells!
        {
            StsException.systemError("StsWellZone.setBot:"
                      + " Zone must have positive thickness.");
        }
        setSubZones();
//        well.setWellLineNeedsRebuild();
    }
    public StsSurfaceVertex getBot(){ return botVertex; }
    public float getBotZ()
    {
        if (botVertex==null) StsException.systemError("StsWellZone.getBotZ:"
                          + " Null top vertex.");
        return botVertex.getPoint().getZ();
    }

    /** set/get flag for display of associated well zones on well lines */
    public void setDisplayOnWellLine(boolean state, boolean setParent)
    {
        displayOnWellLine = state;
        if (wellZoneSet!=null && setParent) wellZoneSet.setDisplayOnWellLines(state);
    }
    public boolean getDisplayOnWellLine() { return displayOnWellLine; }

    public void setWellZoneSet(StsWellZoneSet wellZoneSet)
    {
        if (wellZoneSet==null) StsException.systemError( "StsWellZone.wellZoneSet:  " +
                                                         " well zone set cannot be null.");
        this.wellZoneSet = wellZoneSet;
    }

    public StsWellZoneSet getWellZoneSet() { return wellZoneSet; }

    /** convenience methods */
    public void setName(String name) { wellZoneSet.setName(name); }
    public String getName() { return wellZoneSet.getName(); }
    public boolean setNSubZones(int nSubZones, boolean setParent)
    {
        return wellZoneSet.setNSubZones(nSubZones, setParent);
    }
    public int getNSubZones() { return wellZoneSet.getNSubZones(); }
    public boolean setSubZoneType(int type, boolean setParent)
    {
        return wellZoneSet.setSubZoneType(type, setParent);
    }
    public int getSubZoneType() { return wellZoneSet.getSubZoneType(); }
    public void setStsColor(StsColor color, boolean setParent)
    {
        wellZoneSet.setStsColor(color, setParent);
    }
    public StsColor getStsColor() { return wellZoneSet.getStsColor(); }

    public boolean setZoneType(int type) { return wellZoneSet.setZoneType(type); }
    public int getZoneType() { return wellZoneSet.getZoneType(); }

    // set subzone vector (top to bottom)
    private boolean setSubZones()
    {
        float zoneTop, zoneBot;
        try
        {
            zoneTop = getTopZ();
            zoneBot = getBotZ();
            if(zoneTop >= zoneBot) return false;  // invalid thickness
            rebuildSubZoneZs(zoneTop, zoneBot);
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsWellZone.setSubZones() failed.",
                e, StsException.WARNING);
            return false;
        }
    }

    /* build sub-zone z array, either initially or if nSubZones changed */
    private void rebuildSubZoneZs(float zoneTop, float zoneBot)
    {
        int nSubZones = getNSubZones();
        if (nSubZones == 0)
        {
            subZoneZs = null;
            return;
        }
        int nSubZoneZs = (subZoneZs==null) ? 0 : subZoneZs.length;
        if (nSubZones+1 == nSubZoneZs) return;
//        if (getSubZoneType()==StsZone.SUBZONE_UNIFORM)
        {
            subZoneZs = new float[nSubZones+1];
            float dZ = (zoneBot-zoneTop) / (nSubZones-1);
            for (int i=0; i<nSubZones+1; i++)
            {
                subZoneZs[i] = zoneTop + dZ*i;
            }
        }
    }

    /** return array of subzone boundaries (size is nSubZones+1) */
    public float[] getSubZoneZs()
    {
        if (!setSubZones()) return null;
        return subZoneZs;
    }

    /** return array of layer thicknesses */
    public float[] getLayerThicknesses()
    {
        return null;
    }

    /** return array of subzone thicknesses */
    public float[] getSubZoneThicknesses()
    {
        int nSubZones = getNSubZones();
        if (nSubZones==0 || getSubZoneZs()==null) return null;
        float[] thicknesses = new float[nSubZones];
        for (int i=0; i<nSubZones; i++)
        {
            thicknesses[i] = subZoneZs[i+1]-subZoneZs[i];
        }
        return thicknesses;
    }

    /** return array of subzone center Z values */
    public float[] getLayerCenterZs()
    {
        return null;
    }

    /** return array of subzone center Z values */
    public float[] getSubZoneCenterZs()
    {
        int nSubZones = getNSubZones();
        if (nSubZones==0 || getSubZoneZs()==null) return null;
        float[] centerZs = new float[nSubZones];
        for (int i=0; i<nSubZones; i++)
        {
            centerZs[i] = (subZoneZs[i+1]+subZoneZs[i])/2.0f;
        }
        return centerZs;
    }

    /** get vector of subzone averages for a particular log curve */
    public float[] getLayerAverages(String logCurveName)
    {
        return null;
    }
    public float[] getLayerAverages(StsLogCurve curve)
    {
        return null;
    }

    /** get vector of subzone averages for a particular log curve */
    public float[] getSubZoneAverages(String logCurveName)
    {
        return getSubZoneAverages(well.getLastLogCurveOfType(logCurveName));
    }
    public float[] getSubZoneAverages(StsLogCurve curve)
    {
        // see if we have the curve
        if (curve==null) return null;
        StsObjectRefList logCurves = well.getLogCurves();
        if (logCurves==null) return null;
        if (!logCurves.contains(curve)) return null;

        // get the subzone averages
        int nSubZones = getNSubZones();
        if (nSubZones==0 || getSubZoneZs()==null) return null;
        float[] subZoneAverages = new float[nSubZones];
        for (int i=0; i<nSubZones; i++)
        {
            subZoneAverages[i] = curve.getAverageOverZRange(subZoneZs[i], subZoneZs[i+1]);
        }
        return subZoneAverages;
    }

    /** get array of most commonly occurring category by subzone for a log curve */
    public float[] getSubZoneCategoricalAverages(StsLogCurve curve,
            StsCategoricalFacies categoricalFacies)
    {
        // verify the log curve
        try { if (!well.getLogCurves().contains(curve)) return null; }
        catch (NullPointerException e) { return null; }

        // get the categorical mode values
        int nSubZones = getNSubZones();
        if (nSubZones==0 || getSubZoneZs()==null) return null;
        float[] subZoneAverages = new float[nSubZones];
        for (int i=0; i<nSubZones; i++)
        {
            subZoneAverages[i] = curve.getCategoricalValueOverZRange(categoricalFacies,
                    subZoneZs[i], subZoneZs[i+1]);
        }
        return subZoneAverages;
    }

    /** get top z value to use in ordering well zones */
    public float getOrderingValue()
    {
        try { return topVertex.getPoint().getZ(); }
        catch (Exception e) { return StsParameters.nullValue; }
    }
/*
    private void computeZoneSectionPoints()
    {
        StsGridSectionPoint topPoint, botPoint, edgePoint;

        try
        {
            boolean persistent = true;
            topPoint = topVertex.getSurfacePoint();
            botPoint = botVertex.getSurfacePoint();

            StsSection[] associatedSections = well.getAssociatedSections();
            if(associatedSections == null || associatedSections.length == 0)
            {
                StsException.systemError("StsWellZone.computeZoneSectionPoints() failed." +
                    " No sections connected with this well: " + well.getLabel());
                return;
            }

            // sections[0] will be the section the well is on if not null; otherwise a connected section
            edgePoints = associatedSections[0].constructWellEdgePoints(this, well, topPoint, botPoint);
            int nEdgePoints = edgePoints.getSize();
            // for the remaining connectedSections, add the sectionRowCols for each except for first and last points
            // which should already have them
            for(int n = 1; n < associatedSections.length; n++)
            {
                StsSection associatedSection = (StsSection)associatedSections[n];
                for(int p = 1; p < nEdgePoints-1; p++)
                {
                    edgePoint = (StsGridSectionPoint)edgePoints.getElement(p);
                    edgePoint.addSectionRowCol(associatedSection);
                }
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsSection.sectionPointsOnWell() failed.",
                e, StsException.WARNING);
        }
    }
*/
    /*
    private void constructEdgePoints()
    {
        StsGridSectionPoint topSectionPoint, botSectionPoint;
        StsList sectionPoints;
        StsGridSectionPoint prevSectionPoint, sectionPoint;

        try
        {
            boolean persistent = true;

            topSectionPoint = topVertex.getSurfacePoint();
            botSectionPoint = botVertex.getSurfacePoint();

            float zMin = topSectionPoint.getXYZorT()[2];
            float zMax = botSectionPoint.getXYZorT()[2];

            StsProject project = currentModel.getProject();

            int rowMin = project.getIndexBelow(zMin);
            int rowMax = project.getIndexAbove(zMax);

            int nRows = 0;
            if(rowMax >= rowMin) nRows = rowMax - rowMin + 1;

            sectionPoints = new StsList(nRows + 2, 10);

            // add top and bot points and any section row crossings
            sectionPoints.add(topSectionPoint);

			if(!well.areAllSectionsVertical())
			{
				for(int row = rowMin; row <= rowMax; row++)
				{
					float z = project.getZAtIndex(row);
					StsPoint wellPoint = well.getPointAtZ(z, true);
					sectionPoint = new StsGridSectionPoint(wellPoint, well);
					sectionPoint.setSectionRow();
					sectionPoints.add(sectionPoint);
				}
			}
            sectionPoints.add(botSectionPoint);


            // add section col crossings for section this well is on and sectionRowCols for connected sections
            sectionPoints = addSectionColCrossings(sectionPoints,  persistent);

			// Add points on grid crssoings
			sectionPoints = addGridCrossings(sectionPoints, persistent);

            edgePoints = sectionPoints.convertListToRefList(currentModel, "edgePoints", this);
        }
        catch(Exception e)
        {
            StsException.outputException("StsWellZone.constructEdgePoints() failed.",
                e, StsException.WARNING);
        }
    }
*/
/*
    private void addSectionRowCrossings(StsGridSectionPoint sectionPoint, int row, StsSection[] sections)
    {
        try
        {
            if(sections == null) return;
            int nSections = sections.length;
            for(int n = 0; n < nSections; n++)
            {
                StsSectionRowCol sectionRowCol = sectionPoint.addSectionRowCol(sections[n]);
                if(!sections[n].isVertical()) sectionRowCol.setRowF((float)row);
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsWellZone.addSectionRowCrossings() failed.",
                e, StsException.WARNING);
        }
    }
*/
    private StsList addSectionColCrossings(StsList sectionPoints, boolean persistent)
    {
        StsGridSectionPoint prevSectionPoint, sectionPoint;

        try
        {
            StsSection onSection = well.getOnSection();
            if(onSection != null)
			{
				int nSectionPoints = sectionPoints.getSize();
				StsList newSectionPoints = new StsObjectList(nSectionPoints, 10);

				sectionPoint = (StsGridSectionPoint)sectionPoints.getElement(0);
				newSectionPoints.add(sectionPoint);

				for(int n = 1; n < nSectionPoints; n++)
				{
					prevSectionPoint = sectionPoint;
					sectionPoint = (StsGridSectionPoint)sectionPoints.getElement(n);
					StsList gridPoints = getColCrossings(prevSectionPoint, sectionPoint, onSection, persistent);
					newSectionPoints.addList(gridPoints);
					newSectionPoints.add(sectionPoint);
				}
				sectionPoints = newSectionPoints;
			}

			StsSection[] connectedSections = well.getAllSections();
			for(int n = 0; n < connectedSections.length; n++)
			{
				StsSection section = connectedSections[n];
				int col = 1;
//				int col = section.getWellIndex(well);
				if(col < 0) continue;
				int nSectionPoints = sectionPoints.getSize();
				for(int i = 0; i < nSectionPoints; i++)
				{
					sectionPoint = (StsGridSectionPoint)sectionPoints.getElement(i);
					sectionPoint.setRowOrColIndex(section, StsParameters.COL, col, false);
				}
			}
            return sectionPoints;
        }
        catch(Exception e)
        {
            StsException.outputException("StsWellZone.addSectionRowColCrossings() failed.",
                e, StsException.WARNING);
            return sectionPoints;
        }
    }

    private StsList getColCrossings(StsGridSectionPoint point0, StsGridSectionPoint point1,
                                    StsSection section, boolean persistent)
    {
        float colF0, colF1;
        int col0, col1, inc, nCols, col, n;
        StsGridSectionPoint point;

        colF0 = point0.getSectionRowCol(section).getColF();
        colF1 = point1.getSectionRowCol(section).getColF();

        if(colF1 > colF0)
        {
            col0 = StsMath.above(colF0);
            col1 = StsMath.below(colF1);
            inc = 1;
            nCols = col1 - col0 + 1;
        }
        else
        {
            col0 = StsMath.below(colF0);
            col1 = StsMath.above(colF1);
            inc = -1;
            nCols = col0 - col1 + 1;
        }

        if(nCols <= 0) return null;

        StsList colPoints = new StsList(nCols);
        for(n = 0, col = col0; n < nCols; n++, col += inc)
        {
            double f = (double)((col - colF0)/(colF1 - colF0));
            point = StsGridSectionPoint.sectionInterpolate(point0, point1, f, this, section, persistent);
			point.setRowOrColIndex(section, StsParameters.COL, col);
//            point.getSectionRowCol(section).setColF((float)col);
            colPoints.add(point);
        }
        return colPoints;
    }


    private StsList addGridCrossings(StsList sectionPoints, boolean persistent)
    {
        StsGridSectionPoint sectionPoint, prevSectionPoint;
        StsList gridPoints;

        try
        {
            if(sectionPoints == null) return null;
			if(!well.hasConnectedSection()) return sectionPoints;

            int nSectionPoints = sectionPoints.getSize();
            StsObjectList newSectionPoints = new StsObjectList(nSectionPoints, 10);

            sectionPoint = (StsGridSectionPoint)sectionPoints.getElement(0);
            newSectionPoints.add(sectionPoint);

			StsSection onSection = well.getOnSection();
			StsSection[] connectedSections = well.getAllSections();
			int nConnectedSections = connectedSections.length;

            for(int n = 1; n < nSectionPoints; n++)
            {
                prevSectionPoint = sectionPoint;
                sectionPoint = (StsGridSectionPoint)sectionPoints.getElement(n);

                gridPoints = StsGridSectionPoint.getGridCrossings(prevSectionPoint, sectionPoint, onSection, persistent, NONE);

				if(gridPoints.getSize() > 0)
				{
					for(int s = 0; s < nConnectedSections; s++)
					{
						StsSection connectedSection = connectedSections[s];
						addConnectedSectionRowCols(gridPoints, connectedSection);
					}
					newSectionPoints.addList(gridPoints);
				}
                newSectionPoints.add(sectionPoint);
            }
            return newSectionPoints;
        }
        catch(Exception e)
        {
            StsException.outputException("StsWellZone.addGridCrossings() failed.",
                e, StsException.WARNING);
            return sectionPoints;
        }
    }

    private void addConnectedSectionRowCols(StsList gridPoints, StsSection section)
    {
        StsGridSectionPoint gridPoint;

        if(section == null || gridPoints == null) return;

        try
        {
//			int col;
//			if(well == section.getFirstWell())
//				col = 0;
//			else
//				col = section.getColMax();

            int nGridPoints = gridPoints.getSize();
            for(int n = 0; n < nGridPoints; n++)
            {
                gridPoint = (StsGridSectionPoint)gridPoints.getElement(n);
//                gridPoint.addSectionRowCol(section, this);
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsWellZone.addSectionRowCols() failed.",
                e, StsException.WARNING);
        }
    }

    public StsList getGridEdgePointsList()
    {
        return getEdgePointsList();
    }

    public StsList getEdgePointsList()
    {
        return edgePoints.getList();
    }

    public StsGridSectionPoint getFirstGridPoint()
    {
        if(edgePoints == null) return null;
        else return (StsGridSectionPoint)edgePoints.getFirst();
    }

    public StsGridSectionPoint getLastGridPoint()
    {
        if(edgePoints == null) return null;
        else return (StsGridSectionPoint)edgePoints.getLast();
    }

    public float[][] getXYZPoints()
    {
        StsList points = getEdgePointsList();
        if(points == null) return null;
        int nPoints = points.getSize();
        float[][] xyzPoints = new float[nPoints][];
        for(int n = 0; n < nPoints; n++)
        {
            StsGridSectionPoint gridPoint = (StsGridSectionPoint)points.getElement(n);
            xyzPoints[n] = gridPoint.getPoint().getPointValues();
        }
        return xyzPoints;
    }

    public boolean delete()
    {
        topVertex.deleteEdgeFromVertex(this);
        botVertex.deleteEdgeFromVertex(this);
        if(edgePoints != null) edgePoints.deleteAll();
        super.delete();
        return true;
    }

    public String getLabel()
    {
        return new String("wellEdge-" + getIndex() + " on well " + well.getLabel());
    }
}














