package com.Sts.Actions.Wizards.SimulationFile;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Rescue.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.util.*;

public class StsExportRescueFile
{
    String filename;
    StsModel model;
    StsStatusPanel status;

    static public final boolean rescueOutputBinary = true;

    StsGridDefinition grid;
    int rowMin, rowMax, colMin, colMax; // replace these with a gridRectangle (transient!)
    int cellRowMin, cellRowMax, cellColMin, cellColMax;
    int nCellRows, nCellCols;
    double xMin, xMax, yMin, yMax, xInc, yInc;
    double zMin, zMax;
    int nTotalLayers;

    double kxMultiplier = 1.0;
    double kyMultiplier = 1.0;
    double kzMultiplier = .01;

    RescueCoordinateSystem penUltimateCs = null;
    StsZoneClass zoneClass;  // zones are ordered from the top down

    ArrayList rescueBlocks = new ArrayList();
    Hashtable rescueHorizons = new Hashtable();
    Hashtable rescueSections = new Hashtable();

    StsZoneBlock[][] zoneBlockArray;
    int[][][] actNumArray;

    static final float nullZValue = -99999.0f;


    public StsExportRescueFile(StsModel model, String name, StsStatusPanel status)
    {
        this.status = status;
        this.model = model;
        filename = name;

        zoneClass = (StsZoneClass)model.getCreateStsClass(StsZone.class);
        StsZone topZone = (StsZone)zoneClass.getFirst();
        StsModelSurface topSurface = topZone.getTopModelSurface();
        StsEdgeLoopRadialLinkGrid linkedGrid = topSurface.getLinkedGrid();

        // these are grid (not cell) row-col ranges
        rowMin = linkedGrid.getRowMin();
        rowMax = linkedGrid.getRowMax();
        colMin = linkedGrid.getColMin();
        colMax = linkedGrid.getColMax();

        cellRowMin = rowMin;
        cellRowMax = rowMax-1;
        cellColMin = colMin;
        cellColMax = colMax-1;

        // these are number of cell rows and columns
        nCellRows = rowMax - rowMin;
        nCellCols = colMax - colMin;

        nTotalLayers = zoneClass.getNLayers();

        grid = model.getGridDefinition();

        /** @todo need to generate absolute coordinates in rotated grid */
        double projectXMin = 0.0;
        double projectYMin = 0.0;

//        double projectXMin = model.getProject().getAbsoluteXMin();
//        double projectYMin = model.getProject().getAbsoluteYMin();

        xMin = grid.getXMin() + projectXMin;
        xMax = xMin + grid.getXSize();
        yMin = grid.getYMin() + projectYMin;
        yMax = yMin + grid.getYSize();
        xInc = grid.getXInc();
        yInc = grid.getYInc();
        zMin = model.getProject().getZorTMin();
        zMax = model.getProject().getZorTMax();
    }

    private RescueModel createModel()
    {
        RescueCoordinateSystem ultimateCs = new
            RescueCoordinateSystem("Austrailian TM 146 E",
            RescueCoordinateSystem.LDF, null,
            "northing", "m", "easting", "m", "elevation", "m");

        RescueVertex modelVertex = new RescueVertex("model vertex", ultimateCs, 96.5, 27.2, 1000);

        penUltimateCs = new RescueCoordinateSystem("penultimate cs",
              RescueCoordinateSystem.LDF, modelVertex,
              "northing", "m", "easting", "m", "elevation", "m");

        RescueModel rescueModel = new RescueModel(filename, penUltimateCs);

        return rescueModel;
    }

    private void createBlocks(RescueModel rescueModel)
    {
        StsZone topZone = (StsZone) zoneClass.getFirst();
        StsObject[] zoneBlocks = topZone.getZoneBlocks();
        for( int i=0; i<zoneBlocks.length; i++ )
        {
            StsZoneBlock zb = (StsZoneBlock) zoneBlocks[i];
            String name = "Block" + zb.getName() + i;
            RescueBlock block = new RescueBlock(name, rescueModel);
            rescueBlocks.add(block);
        }
    }

    private void createHorizons(RescueModel rescueModel)
    {
        StsMainObject[] horizons = model.getObjectListOfType(StsModelSurface.class, StsModelSurface.MODEL);
        for( int i=0; i<horizons.length; i++ )
        {
            String name = "Horizon-" + horizons[i].getName();
            RescueHorizon h = new RescueHorizon(name, rescueModel);
            rescueHorizons.put(horizons[i], h);
        }
    }

    private void createSections(RescueModel rescueModel)
    {
        StsSection[] sections = (StsSection[])model.getCastObjectList(StsSection.class);
        for( int i=0; i<sections.length; i++ )
        {
            StsSection section = sections[i];
            int type = 0;
            switch( section.getType() )
            {
                case StsSection.AUXILIARY:
                    type = RescueSection.AUXILLIARY;
                    break;
                case StsSection.BOUNDARY:
                    type = RescueSection.LEASE_BOUNDARY;
                    break;
                case StsSection.FAULT:
                    type = RescueSection.FAULT;
            }
            RescueSection rescueSection = new RescueSection(RescueCoordinateSystem.LDF,
                section.getName(), rescueModel, type,
                0, section.getNCols(),
                0, section.getNRows(),
                StsParameters.nullValue);

            rescueSections.put(section, section);

            rescueSection.Geometry().AssignXValue(section.getPatch().getXValues());
            rescueSection.Geometry().AssignYValue(section.getPatch().getYValues());
            rescueSection.Geometry().AssignZValue(section.getPatch().getZValues());
        }
    }

    private void createUnits(RescueModel rescueModel)
    {
        ArrayList units = new ArrayList();
        StsObject[] zones = zoneClass.getElements();
        for( int i=0; i<zones.length; i++ )
        {
            StsZone zone = (StsZone) zones[i];
            RescueUnit unit = new RescueUnit(zone.getName(), rescueModel);
            units.add(unit);

            RescueHorizon top = (RescueHorizon) rescueHorizons.get(zone.getTopModelSurface());
            RescueHorizon base = (RescueHorizon) rescueHorizons.get(zone.getBaseModelSurface());
            unit.SetHorizonAboveMe(top);
            unit.SetHorizonBelowMe(base);
            top.SetUnitBelowMe(unit);
            base.SetUnitAboveMe(unit);

            // create all of the block units for this zone
            createBlockUnits(zone, rescueModel, unit, rescueBlocks);
        }
    }


    private void createBlockUnits(StsZone zone, RescueModel rescueModel, RescueUnit unit, List blocks)
    {
        StsObject[] zoneBlocks = zone.getZoneBlocks();
        for( int i=0; i<zoneBlocks.length; i++ )
        {
            StsZoneBlock zb = (StsZoneBlock) zoneBlocks[i];
            RescueBlock block = (RescueBlock) rescueBlocks.get(i);
            RescueBlockUnit blockUnit = createBlockUnit(zb, rescueModel, block, unit);
        }
    }

    private RescueBlockUnit createBlockUnit(StsZoneBlock zoneBlock,
                                            RescueModel rescueModel,
                                            RescueBlock block,
                                            RescueUnit unit)
    {
        zoneBlock.constructLayerGrids();

        StsZone zone = zoneBlock.getZone();
        int nCols = zoneBlock.getNCols();
        int nRows = zoneBlock.getNRows();
        int nz = zone.getNLayerGrids();

        RescueVertex offset = new RescueVertex("model offset", penUltimateCs, 0,0,0);
        RescueBlockUnit myBlockUnit = new RescueBlockUnit(RescueCoordinateSystem.LDF,
                block, unit,
                zoneBlock.xMin, zoneBlock.xInc, zoneBlock.colMin, zoneBlock.colMax,
                zoneBlock.yMin, zoneBlock.yInc, zoneBlock.rowMin, zoneBlock.rowMax,
                0, nz,
                StsParameters.nullValue,  //missing value
                0.0f, offset);

        // Set Z Value Geometry for Block Unit
        myBlockUnit.GridGeometry().AssignZValue(getBlockUnitGeometry(zoneBlock, nCols, nRows, nz));


        // add the bounding surfaces
        RescueHorizon horizon = (RescueHorizon) rescueHorizons.get(zone.getBaseModelSurface());
        RescueBlockUnitHorizonSurface bottom = createBlockUnitSurface(horizon, zoneBlock.getBottomGrid());
        myBlockUnit.SetSurfaceBelowMe(bottom);
        bottom.SetBlockUnitAboveMe(myBlockUnit);

        horizon = (RescueHorizon) rescueHorizons.get(zone.getTopModelSurface());
        RescueBlockUnitHorizonSurface top = createBlockUnitSurface(horizon, zoneBlock.getTopGrid());
        myBlockUnit.SetSurfaceAboveMe(top);
        top.SetBlockUnitBelowMe(myBlockUnit);

        // add the boundary
        Hashtable sides = createBlockUnitSides(rescueModel, zoneBlock);
        Enumeration elements = sides.elements();
        while( elements.hasMoreElements() )
        {
            RescueBlockUnitSide side = (RescueBlockUnitSide) elements.nextElement();
            myBlockUnit.MacroVolume().AddBlockUnitSide(side);
        }

        // add the k trim loops
        createEdgeSets(rescueModel, myBlockUnit, top, bottom, zoneBlock);

        // add the face loops
        createFaceLoops(rescueModel, myBlockUnit, top, bottom, sides, zoneBlock);

        return myBlockUnit;
    }

    private Hashtable createBlockUnitSides(RescueModel rescueModel, StsZoneBlock zoneBlock)
    {
        Object[] zoneSides = zoneBlock.getSides();
        Hashtable sides = new Hashtable(zoneSides.length);
        for( int i=0; i<zoneSides.length; i++ )
        {
            StsZoneSide zoneSide = (StsZoneSide) zoneSides[i];
            RescueSection rescueSection = (RescueSection) rescueSections.get(zoneSide.getSection());
            RescueBlockUnitSide side = new RescueBlockUnitSide(rescueSection);
            sides.put(zoneSide, side);

//            section.Geometry().AssignXValue(zoneSide.getXValues());
//            section.Geometry().AssignYValue(zoneSide.getYValues());
//            section.Geometry().AssignZValue(zoneSide.getZValues());
        }
        return sides;
    }

    private RescueBlockUnitHorizonSurface createBlockUnitSurface(RescueHorizon rescueHorizon, StsBlockGrid grid)
    {
        int nCols = grid.getNCols();
        int nRows = grid.getNRows();
        RescueBlockUnitHorizonSurface surface = new RescueBlockUnitHorizonSurface(
            RescueCoordinateSystem.LDF, rescueHorizon,
            grid.getXMin(), (float)xInc, grid.getColMin(), nCols,
            grid.getYMin(), (float)yInc, grid.getRowMin(), nRows,
//            grid.getXMin(), grid.getXInc(), grid.getColMin(), nCols,
//            grid.getYMin(), grid.getYInc(), grid.getRowMin(), nRows,
            StsParameters.nullValue, grid.getAngle(), null, RescueSurface.HORIZON);

        surface.Geometry().AssignZValue(getBlockUnitHorizonSurfaceGeometry(grid, nCols, nRows));
        return surface;
    }

    private float[] getBlockUnitHorizonSurfaceGeometry(StsBlockGrid grid, int nCols, int nRows)
    {
        float[] surfaceGeometry = new float[nCols*nRows];
        int n = 0;
        for( int i=0; i<nRows; i++ )
            for( int j=0; j<nCols; j++, n++ )
                surfaceGeometry[n] = grid.getPointZ(i, j);
        return surfaceGeometry;
    }

    private float[] getBlockUnitGeometry(StsZoneBlock zoneBlock,
                                            int nCols, int nRows, int nz)
    {
        float[] blockUnitGeometry = new float[nCols*nRows*nz];
        int n = 0;
        for( int k=0; k<nz; k++ )
            for( int i=0; i<nRows; i++ )
                for( int j=0; j<nCols; j++, n++ )
                    blockUnitGeometry[n] = zoneBlock.getSubZoneGridZ(k, i, j);
        return blockUnitGeometry;
    }


    private List createEdgeSets(RescueModel rescueModel, RescueBlockUnit bu,
                                RescueBlockUnitHorizonSurface top,
                                RescueBlockUnitHorizonSurface bottom,
                                StsZoneBlock zb)
    {
        ArrayList sets = new ArrayList();
        int nLayers = zb.getZone().getNSubZones() + 1;
        for( int k=0; k<nLayers; k++ )
        {
            RescueEdgeSet set = createEdgeSet(rescueModel, k, bu, top, bottom, zb);
            bu.MacroVolume().AddKLayerEdge(set); // Is this redundant for 0 and nzcorn -1  Apparently not

            if( k == 0 ) bu.MacroVolume().SetTopEdge(set);
            else if( k == nLayers-1 ) bu.MacroVolume().SetBottomEdge(set);

            sets.add(set);
        }
        return sets;
    }

    private RescueEdgeSet createEdgeSet(RescueModel rescueModel, int layer,
                                        RescueBlockUnit bu,
                                        RescueBlockUnitHorizonSurface top,
                                        RescueBlockUnitHorizonSurface bottom,
                                        StsZoneBlock zb)
    {
        int nLayers = zb.getZone().getNSubZones() + 1;
        RescueEdgeSet set = new RescueEdgeSet();
        RescueTrimLoop trimLoop = new RescueTrimLoop();

        // build a lookup for section side via edge
        Object[] sideArray = zb.getSides();
        Hashtable edgeToSectionMap = new Hashtable(sideArray.length);
        for( int i=0; i<sideArray.length; i++ )
        {
            StsZoneSide side = (StsZoneSide) sideArray[0];
            edgeToSectionMap.put(side.getTopEdge(), side.getSection());
        }

        StsEdgeLoop loop = zb.getTopGrid().getEdgeLoop();
        Object[] edges = loop.getDirectedEdges().getList();
        for( int i=0; i<edges.length; i++ )
        {
            StsDirectedEdge edge = (StsDirectedEdge) edges[i];

            // first point
            StsGridSectionPoint firstPoint = edge.getFirstPoint();
            float[] xyz = firstPoint.getXYZorT();
//            float[] xyz = getPointValues(layer, firstPoint, zb);
            RescueTrimVertex firstVertex = new RescueTrimVertex(rescueModel,
                                        xyz[0], xyz[1], xyz[2]);

            // set the UV value on the section
//            StsSection section = firstPoint.getSection();
            StsSection section = (StsSection) edgeToSectionMap.get(edge.getEdge());
            if( section != null )
            {
                RescueSection rescueSection = (RescueSection) rescueSections.get(section);
                firstVertex.SetUVValue(rescueSection, 0f, xyz[2]);
            }

            // set the UV value on the surface
            if( layer == 0 )
                firstVertex.SetUVValue(top, firstPoint.getGridColF(), firstPoint.getGridRowF());
            else if( layer == nLayers-1 )
                firstVertex.SetUVValue(bottom, firstPoint.getGridColF(), firstPoint.getGridRowF());

            // last point
            StsGridSectionPoint lastPoint = edge.getLastPoint();
            xyz = lastPoint.getXYZorT();
//            xyz = getPointValues(layer, lastPoint, zb);
            RescueTrimVertex lastVertex = new RescueTrimVertex(rescueModel,
                                        xyz[0], xyz[1], xyz[2]);

            // set the UV value on the section
//            section = lastPoint.getSection();
            if( section != null )
            {
                RescueSection rescueSection = (RescueSection) rescueSections.get(section);
                lastVertex.SetUVValue(rescueSection, 1f, xyz[2]);
            }

            // set the UV value on the surface
            if( layer == 0 )
                lastVertex.SetUVValue(top, lastPoint.getGridColF(), lastPoint.getGridRowF());
            else if( layer == nLayers-1 )
                lastVertex.SetUVValue(bottom, lastPoint.getGridColF(), lastPoint.getGridRowF());

            RescuePolyLine polyLine = new RescuePolyLine(rescueModel, firstVertex, lastVertex);

            // now add all the points in between the first and last
            StsList edgePoints = edge.getEdgePoints();
            for( int j=1; j<edgePoints.getSize()-1; j++ )
            {
                StsGridSectionPoint p = (StsGridSectionPoint) edgePoints.getElement(j);
                xyz = p.getXYZorT();
//                xyz = getPointValues(layer, p, zb);
                RescuePolyLineNode node = new RescuePolyLineNode(xyz[0], xyz[1], xyz[2]);
//                section = p.getSection();
                if( section != null )
                {
                    RescueSection rescueSection = (RescueSection) rescueSections.get(section);
                    node.SetUVValue(rescueSection, j / edgePoints.getSize(), xyz[2]);
                }

                if( layer == 0 )
                    node.SetUVValue(top, p.getGridColF(), p.getGridRowF());
                else if( layer == nLayers-1 )
                    node.SetUVValue(bottom, p.getGridColF(), p.getGridRowF());

                polyLine.AddPolyLineNode(node);
            }

            // determine the direction
            int direction = RescueTrimEdge.R_LEFT_TO_RIGHT;
            switch( edge.getDirection() )
            {
            case StsDirectedEdge.MINUS:
                direction = RescueTrimEdge.R_LEFT_TO_RIGHT;
            case StsDirectedEdge.PLUS:
                direction = RescueTrimEdge.R_RIGHT_TO_LEFT;
            }

            // finally, add the edge to the loop
            trimLoop.AddLoopEdge(new RescueTrimEdge(polyLine, direction));
        }
        set.AddBoundaryLoop(trimLoop);
        return set;
    }
/*  xyz[2] should be the same z that is already in p, so don't need to recompute. TJL 5/26/00
    private float[] getPointValues(int layer, StsGridSectionPoint p, StsZoneBlock zb)
    {
        float[] xyz = p.getPointXYZ();
        if( layer != 0 )
            xyz[2] = zb.getSubZoneGridZ(layer, p.getGridRowF(), p.getGridColF());
        return xyz;
    }
*/

    private List createFaceLoops(RescueModel rescueModel, RescueBlockUnit bu,
                                RescueBlockUnitHorizonSurface top,
                                RescueBlockUnitHorizonSurface bottom,
                                Hashtable sides, StsZoneBlock zb)
    {
        Hashtable surfaces = new Hashtable(2);
        surfaces.put(zb.getZone().getBaseModelSurface(), bottom);
        surfaces.put(zb.getZone().getTopModelSurface(), top);

        Enumeration keys = sides.keys();
        while( keys.hasMoreElements() )
        {
            StsZoneSide zoneSide = (StsZoneSide) keys.nextElement();
            RescueTrimLoop loop = createFaceLoop(rescueModel, bu, surfaces, zb, zoneSide);

            RescueBlockUnitSide side = (RescueBlockUnitSide) sides.get(zoneSide);
            side.Edges().AddBoundaryLoop(loop);
        }

        return null;
    }


    private RescueTrimLoop createFaceLoop(RescueModel rescueModel, RescueBlockUnit bu,
                                            Hashtable surfaces,
                                            StsZoneBlock zb, StsZoneSide side)
    {
        RescueTrimLoop loop = new RescueTrimLoop();

        StsEdgeLoop edgeLoop = side.getEdgeLoop();

        Enumeration edges = edgeLoop.getDirectedEdges().elements();
        while( edges.hasMoreElements() )
        {
            StsDirectedEdge edge = (StsDirectedEdge) edges.nextElement();

            // first point
            StsGridSectionPoint firstPoint = edge.getFirstPoint();
            float[] xyz = firstPoint.getXYZorT();
            RescueTrimVertex firstVertex = new RescueTrimVertex(rescueModel,
                                        xyz[0], xyz[1], xyz[2]);


            // set the UV value on the section
//            StsSection section = firstPoint.getSection();
            StsSection section = side.getSection();
            if( section != null )
            {
                RescueSection rescueSection = (RescueSection) rescueSections.get(section);
                firstVertex.SetUVValue(rescueSection, 0f, xyz[2]);
            }

            // set the UV value on the surface
            StsSurfaceVertex vertex = firstPoint.getVertex();
            if( vertex != null )
            {
                RescueBlockUnitHorizonSurface surface = (RescueBlockUnitHorizonSurface) surfaces.get(vertex.getSurface());
                firstVertex.SetUVValue(surface, firstPoint.getGridColF(), firstPoint.getGridRowF());
            }

            // last point
            StsGridSectionPoint lastPoint = edge.getLastPoint();
            xyz = lastPoint.getXYZorT();
            RescueTrimVertex lastVertex = new RescueTrimVertex(rescueModel,
                                        xyz[0], xyz[1], xyz[2]);

            // set the UV value on the section
//            section = lastPoint.getSection();
            if( section != null )
            {
                RescueSection rescueSection = (RescueSection) rescueSections.get(section);
                lastVertex.SetUVValue(rescueSection, 1f, xyz[2]);
            }

            // set the UV value on the surface
            vertex = lastPoint.getVertex();
            if( vertex != null )
            {
                RescueBlockUnitHorizonSurface surface = (RescueBlockUnitHorizonSurface) surfaces.get(vertex.getSurface());
                lastVertex.SetUVValue(surface, lastPoint.getGridColF(), lastPoint.getGridRowF());
            }

            RescuePolyLine polyLine = new RescuePolyLine(rescueModel, firstVertex, lastVertex);

            // now add all the points in between the first and last
            StsList edgePoints = edge.getEdgePoints();
            for( int j=1; j<edgePoints.getSize()-1; j++ )
            {
                StsGridSectionPoint p = (StsGridSectionPoint) edgePoints.getElement(j);
                xyz = p.getXYZorT();
                RescuePolyLineNode node = new RescuePolyLineNode(xyz[0], xyz[1], xyz[2]);
//                section = p.getSection();
                if( section != null )
                {
                    RescueSection rescueSection = (RescueSection) rescueSections.get(section);
                    node.SetUVValue(rescueSection, j/edgePoints.getSize(), xyz[2]);
                }

                // set the UV value on the surface
                vertex = p.getVertex();
                if( vertex != null )
                {
                    RescueBlockUnitHorizonSurface surface = (RescueBlockUnitHorizonSurface) surfaces.get(vertex.getSurface());
                    node.SetUVValue(surface, p.getGridColF(), p.getGridRowF());
                }

                polyLine.AddPolyLineNode(node);
            }

            // determine the direction
            int direction = RescueTrimEdge.R_LEFT_TO_RIGHT;
            switch( edge.getDirection() )
            {
            case StsDirectedEdge.MINUS:
                direction = RescueTrimEdge.R_LEFT_TO_RIGHT;
            case StsDirectedEdge.PLUS:
                direction = RescueTrimEdge.R_RIGHT_TO_LEFT;
            }

            // finally, add the edge to the loop
            loop.AddLoopEdge(new RescueTrimEdge(polyLine, direction));
        }
        return loop;
    }

    public boolean saveModel(StsStatusUI status)
    {
        try
        {
            status.setMaximum(100);
            System.out.println(RescueModel.BuildDate());
        }
        catch(UnsatisfiedLinkError ule)
        {
            new StsMessage(model.win3d,  StsMessage.WARNING, "Sorry. RESCUE library is not available.");
            return false;
        }

       try
        {
            RescueModel rescueModel = createModel();
            status.setProgress(10);
            createSections(rescueModel);
            status.setProgress(20);
            createBlocks(rescueModel);
            status.setProgress(30);
            createHorizons(rescueModel);
            status.setProgress(40);
            createUnits(rescueModel);
            status.setProgress(80);

            if(rescueOutputBinary)
                rescueModel.ArchiveModel(filename+".bin", rescueOutputBinary);
            else
                rescueModel.ArchiveModel(filename+".txt", rescueOutputBinary);
            status.setProgress(80);



/**************************************************

  Note that in order to avoid holding the entire model in memory at once
  we have to first archive the model, then build the separately archivable
  parts, unload them, then archive again to make sure we got any related
  changes.

******************************************************/

/*
            status.setProgress(10);

            model.UnloadWireframe();
            model.ArchiveModel();
            model.dispose();

            status.setProgress(100);
*/
            return true;
        }
        catch(Exception e)
        {
            status.setProgress(0);
            StsException.outputException("StsExportRescueFile.saveModel() failed.", e, StsException.WARNING);
            return false;
        }
    }

}
