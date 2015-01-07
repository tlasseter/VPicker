
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.MVC.*;
import com.Sts.UI.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;
import com.Sts.Types.*;
import com.Sts.Interfaces.*;

import java.awt.*;
import java.io.*;

public class StsBuiltModel extends StsMainTimeObject implements StsTreeObjectI
{
    private StsObjectRefList zones;
    private StsObjectRefList blocks;
    private String filesPathname; 

    // this is a combined list used by the combo box
    // private transient StsPropertyType displayPropertyType = nonePropertyType;
    transient public StsGridBoundingBox modelBoundingBox;
    transient int nCellLayers;
    transient int nTotalCells;
    transient StsSurfaceEdge[] edgeList;
    transient StsSurfaceEdge pickedSurfaceEdge = null;

    static StsObjectPanel objectPanel = null;
    static public StsFieldBean[] displayFields = null;
    static public StsFieldBean[] propertyFields = null;
    // convenience flag copies
    static final int NOT_BUILDABLE = StsZone.NOT_BUILDABLE;
    static final int BUILDABLE = StsZone.BUILDABLE;
    static final int BUILT = StsZone.BUILT;

    public StsBuiltModel()
    {
    }

   public StsBuiltModel(StsZone[] zoneList, StsBlock[] blockList)
   {
       super(false);
       int nZones = zoneList.length;
       zones = StsObjectRefList.constructor(nZones, 2, "zones", this);
       for(int n = 0; n < nZones; n++)
           zoneList[n].setZoneNumber(n);
       zones.add(zoneList);

       int nBlocks = blockList.length;
       blocks = StsObjectRefList.constructor(nBlocks, 2, "blocks", this);
       blocks.add(blockList);
       addToModel();
       String projectName = currentModel.getProject().getName();
       setName(projectName + "-" + getIndex());
       filesPathname = currentModel.getProject().getRootDirString() + "EclipseFiles" + File.separator + name + File.separator;
    }

    public String getFilesPathname()
    {
        File filesDirectory = new File(filesPathname);
        if(filesDirectory.exists()) return filesPathname;
        filesPathname = currentModel.getProject().getRootDirString() + "EclipseFiles" + File.separator + name + File.separator;
        dbFieldChanged("filesPathname", filesPathname);
        return filesPathname;
    }

    public boolean initialize(StsModel model)
    {
        constructLayers();
        constructBlockLayout();
        return true;
    }

    public void constructLayers()
    {
        StsZone[] zones = getZones();
        nCellLayers = 0;
        for (StsZone zone : zones)
        {
            StsObjectRefList subZones = zone.getSubZones();
            if(subZones == null) continue;
            int nSubZones = subZones.getSize();
            for (int sz = 0; sz < nSubZones; sz++)
            {
                StsSubZone subZone = (StsSubZone) subZones.getElement(sz);
                nCellLayers += subZone.getNLayers();
            }
        }
    }

    private void constructBlockLayout()
    {
        StsBlock[] blocks = getBlocks();
        int nBlocks = blocks.length;
        int nCellRows = 0;
        int nCellCols = 0;
        for (int n = 0; n < nBlocks; n++)
        {
            StsBlock block = blocks[n];
            block.eclipseRowMin = nCellRows + 1;
            nCellRows += block.nCellRows;
            nCellCols = Math.max(nCellCols, block.nCellCols);
            // insert a null row between blocks
            if(n < nBlocks-1) nCellRows++;
            block.setLayerRange(nCellLayers);
        }
        nTotalCells = nCellRows * nCellCols * nCellLayers;
    }

    /** For now, leave it up to the user to build model rather than building it automatically.
	 * This way, user can change parameters before rebuilding.
	 */
   /*
    public boolean initialize(StsModel model)
    {
        return true;
//        StsStatusUI status = StsStatusArea.getStatusArea();
//        return classInitialize(model, false, status);
    }

	public boolean initialize()
	{
        StsStatusPanel statusPanel = StsStatusPanel.constructStatusDialog(currentModel.win3d, "Building model...");
		boolean ok = buildModel(currentModel, true, statusPanel);
        StsStatusPanel.disposeDialog();
        return ok;
	}
    */
 
    public boolean buildModel(StsModel model, StsStatusUI status)
    {
        float progress, maxProgress, incProgress;
        int nZones = zones.getSize();
        if(nZones == 0) return false;
        int nBlocks = blocks.getSize();
        if(nBlocks == 0) return false;

        constructHardWiredSubZones();

        StsCursor cursor = new StsCursor(model.win3d, Cursor.WAIT_CURSOR);

        try
        {
            boolean success;

//            model.disableDisplay();

            status.setMaximum(100.0f);

			model.disableDisplay();

            progress = 0.0f;
            maxProgress = 5.0f;

            status.setProgress(maxProgress);

            for(int n = 0; n < nBlocks; n++)
            {
                StsBlock block = (StsBlock)blocks.getElement(n);
                block.deleteTransientArrays();
            }

            // get surfaces from zones in top to bottom sequence
            StsModelSurfaceClass surfaceClass = (StsModelSurfaceClass)model.getCreateStsClass(StsModelSurface.class);
            StsModelSurface[] surfaces = surfaceClass.getTopDownModelSurfaces();
            int nSurfaces = surfaces.length;
            if(nSurfaces == 0) return false;

            // turn off original surfaces, turn on new surfaces
            for(StsModelSurface surface : surfaces)
            {
                StsSurface origSurface = surface.getOriginalSurface();
                if(origSurface != null && origSurface.checkIsLoaded()) origSurface.unload();
                if(surface.getIsVisible()) surface.setIsVisible(true);
                surface.setType(StsSurface.MODEL);
                surface.deleteTransientArrays();
                surface.initializePointArrays();
            }

            // gap surface grids
            progress = maxProgress;
            maxProgress = 30.0f;
            incProgress = (maxProgress-progress)/nSurfaces;
            success = true;
            for(StsModelSurface surface : surfaces)
            {
                if(!surface.gapSurfaceGrid(this, status, progress, incProgress)) success = false;
                progress += incProgress;
            }
            if(!success) return false;

            // from edgeLoops defined by blockSide intersections with surfaces, compute
            // boundingBoxes for each block
            StsBlock[] blockList = getBlocks();
            modelBoundingBox = StsBlock.computeModelBoundingBox(blockList);
            if(modelBoundingBox == null) return false;

            // rotatedPoints for each sectionLine are currently based on spacing of original definition points
            // reset them so that they are on even-spacing as defined by this builtModel
            resetSectionLinePoints(model, modelBoundingBox);

            // for each surface, construct a blockGrid for each block from the bottom up
            // as we build the next blockGrid above in a block, don't allow it to
            // penetrate the blockGrid below
            progress = maxProgress;
            maxProgress = 90.0f;
            incProgress = (maxProgress-progress)/nSurfaces;
            success = true;
            for(int n = nSurfaces-1; n >= 0; n--)
            {
                if(!surfaces[n].constructBlockGrids(status, progress, incProgress)) success = false;
                progress += incProgress;
            }
            if(!success) return false;

            progress = maxProgress;
            maxProgress = 100.0f;
            incProgress = (maxProgress-progress)/nZones;
            int topLayerNumber = 0;
            // zoneBlocks aren't persistent; reinitialize this colorIndex static so zoneBlocks are colored in consistent sequence
            StsZoneBlock.nextColorIndex = 0;
            for(int n = 0; n < nZones; n++)
            {
                StsZone zone = (StsZone)zones.getElement(n);
                zone.setTopLayerNumber(topLayerNumber);
                if(!zone.buildZoneBlocks(blockList, status, progress, incProgress)) success = false;
                topLayerNumber += zone.getNLayers();
                progress += incProgress;
            }

			// toggle off sectionEdges and toggle on surfaceEdges
            StsSectionClass sectionClass = (StsSectionClass)model.getCreateStsClass(StsSection.class);
			sectionClass.setDisplaySectionEdges(false);
			surfaceClass.setDisplayHorizonEdges(true);

            if(!Main.debugPoint)
            {
                for(int n = 0; n < nSurfaces; n++)
                surfaces[n].deleteTransientArrays();
            }

			model.enableDisplay();
            return true;
//            if(success) model.win3d.win3dDisplay();
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "buildModel", e);
            return false;
        }
        finally
        {
		    // checkCreateTransientBuildModelTransaction(model);
			// if(actionManager != null) model.getDatabase().saveTransactions();
            cursor.restoreCursor();
            model.enableDisplay();
        }
    }

    public boolean rebuildModel(StsModel model, StsStatusPanel statusPanel)
    {
        if(!rebuildSectionsOK()) return false;
        if(!rebuildSurfacesOK()) return false;
        buildModel(currentModel, statusPanel);
        return true;
    }

    private boolean rebuildSectionsOK()
    {
        StsProject project = currentModel.getProject();
        StsObject[] sectionObjects = currentModel.getObjectList(StsSection.class);
        if(project.hasVelocityModel())
        {
            for(int n = 0; n < sectionObjects.length; n++)
            {
                StsSection section = (StsSection)sectionObjects[n];
                section.constructSection();
            }
            return true;
        }
        // if no velocity model, we can rebuild if this is the original zDomain for all sections
        for(int n = 0; n < sectionObjects.length; n++)
        {
            StsSection section = (StsSection)sectionObjects[n];
            if(!section.isZDomainOriginal(project.zDomain)) return false;
        }
        return true;
    }

   private boolean rebuildSurfacesOK()
   {
       StsModelSurfaceClass surfaceClass = (StsModelSurfaceClass)currentModel.getCreateStsClass(StsModelSurface.class);
       StsModelSurface[] surfaces = surfaceClass.getTopDownModelSurfaces();
       if(surfaces.length == 0) return false;
       for(StsModelSurface surface : surfaces)
           surface.initialize(currentModel);
       return true;
   }

    private void constructHardWiredSubZones()
    {
        int nZones = zones.getSize();
        if(nZones == 1)
        {
            int nSubZones = 10;
            for(int n = 0; n < nZones; n++)
            {
                StsZone zone = (StsZone)zones.getElement(n);
                zone.constructHardWiredSubZones(nSubZones);
            }
        }
        else if(nZones == 10)
        {
            int nSubZones = 1;
            for(int n = 0; n < nZones; n++)
            {
                StsZone zone = (StsZone)zones.getElement(n);
                zone.constructHardWiredSubZones(nSubZones);
            }
        }
        else if(nZones == 3)
        {
            StsZone zone;
            zone = (StsZone)zones.getElement(0);
            zone.constructHardWiredSubZones(1);
            zone = (StsZone)zones.getElement(1);
            zone.constructHardWiredSubZones(8);
            zone = (StsZone)zones.getElement(2);
            zone.constructHardWiredSubZones(1);
        }
    }

    private void resetSectionLinePoints(StsModel model, StsGridBoundingBox modelBoundingBox)
    {
        StsLine line;

        try
        {
            StsObject[] sections = model.getObjectList(StsSection.class);
            int nSections = sections.length;
            for(int n = 0; n < nSections; n++)
            {
                StsSection section = (StsSection)sections[n];
                line = section.getFirstLine();
                line.resetSectionLinePoints(modelBoundingBox);

                line = section.getLastLine();
                if(line.onSection != null)
                    line.resetSectionLinePoints(modelBoundingBox);
            }
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "resetSectionLinePoints", e);
        }
    }

    public boolean delete()
    {
		StsActionManager actionManager = null;

        StsCursor cursor = new StsCursor(currentModel.win3d, Cursor.WAIT_CURSOR);
        try
        {
			actionManager = currentModel.mainWindowActionManager;
            // if(actionManager != null) model.getDatabase().blockTransactions();
            // currentModel.deleteStsClass(StsBuiltModel.class);
            // remove all blocks
            // model.deleteStsClass(StsBlock.class);
            // zones remain, but zoneBlocks will be rebuilt
            // StsZoneBlock(s) are not persistent, so we don't need to delete them from the model
            StsZone[] zones = (StsZone[])currentModel.getCastObjectList(StsZone.class);
            for(StsZone zone : zones)
                zone.deleteModel();
            // delete non-persistent surfaces edges: entire class
            currentModel.deleteStsClass(StsSurfaceEdge.class);
            currentModel.deleteNonPersistentObjects(StsGridSectionPoint.class);
            // modelSurfaces contain some model transients; remove them and reset stuff as needed
            StsModelSurfaceClass surfaceClass = (StsModelSurfaceClass)currentModel.getCreateStsClass(StsModelSurface.class);
            StsModelSurface[] modelSurfaces = surfaceClass.getTopDownModelSurfaces();
            int nSurfaces = modelSurfaces.length;
            if(nSurfaces == 0) return false;

            for(int n = 0; n < nSurfaces; n++)
            {
                StsModelSurface modelSurface = modelSurfaces[n];
                if(modelSurface == null) continue;
                modelSurface.deleteModel();
                modelSurface.resetNullTypes();
            }

            StsSectionClass sectionClass = (StsSectionClass)currentModel.getCreateStsClass(StsSection.class);
			sectionClass.setDisplaySectionEdges(true);
            /*
            int nSections = sectionClass.getSize();
            for(int n = 0; n < nSections; n++)
            {
                StsSection section = (StsSection)sectionClass.getElement(n);
                section.constructSection();
            }
            */
            super.delete();
            currentModel.refreshObjectPanel();
            return true;
        }
        catch(Exception e)
        {
            StsException.outputWarningException(StsBuiltModel.class, "deleteModel", e);
            return false;
        }
        finally
        {
             cursor.restoreCursor();
			 if(actionManager != null) currentModel.getDatabase().saveTransactions();
             currentModel.enableDisplay();
        }
    }

    static public void checkDeleteBoundary(StsModel model)
    {
        model.deleteObjectsOfType(StsSection.class, StsParameters.BOUNDARY);
 //       model.deleteObjectsOfType(StsLine.class, StsParameters.BOUNDARY);

    }

    static public void checkDeleteFaults(StsModel model)
    {
        model.deleteObjectsOfType(StsSection.class, StsParameters.FAULT);
    }

    public StsZone[] getZones()
    {
        return (StsZone[])zones.getCastList();
    }

    public StsBlock[] getBlocks()
    {
        return (StsBlock[])blocks.getCastList();
    }

    public StsSurfaceEdge[] buildSurfaceEdgeList()
    {
        if(blocks == null) return null;
        edgeList = new StsSurfaceEdge[100];
        int nEdges = 0;
        int nBlocks = blocks.getSize();
        for(int b = 0; b < nBlocks; b++)
        {
            StsBlock block = (StsBlock)blocks.getElement(b);
            StsList blockGrids = block.blockGrids;
            int nBlockGrids = blockGrids.getSize();
            for(int g = 0; g < nBlockGrids; g++)
            {
                StsBlockGrid blockGrid = (StsBlockGrid)blockGrids.getElement(g);
                StsList edges = blockGrid.getEdges();
                int nGridEdges = edges.getSize();
                for(int e = 0; e < nGridEdges; e++)
                {
                    StsSurfaceEdge edge = (StsSurfaceEdge)edges.getElement(e);
                    edge.setIndex(nEdges);
                    edgeList = (StsSurfaceEdge[])StsMath.arrayAddElement(edgeList, edge, nEdges, 100);
                    nEdges++;
                }
            }
        }
        edgeList = (StsSurfaceEdge[])StsMath.trimArray(edgeList, nEdges);
        return edgeList;
    }

    public void setPickedSurfaceEdge(StsSurfaceEdge surfaceEdge)
    {
        pickedSurfaceEdge = surfaceEdge;
    }

    public boolean anyDependencies()
    {
        return false;
    }

    public StsFieldBean[] getDisplayFields()
    {
        return null;
    }

    public StsFieldBean[] getPropertyFields()
    {
        return propertyFields;
    }

    public Object[] getChildren()
    {
        return new Object[0];
    }

    public StsObjectPanel getObjectPanel()
    {
        if (objectPanel == null) objectPanel = StsObjectPanel.constructor(this, true);
        return objectPanel;
    }

    public void treeObjectSelected()
    {
        currentModel.getCreateStsClass("com.Sts.DBTypes.StsBuiltModel").selected(this);
    }
}
