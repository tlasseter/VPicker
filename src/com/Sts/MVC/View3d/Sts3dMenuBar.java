
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System


package com.Sts.MVC.View3d;

import com.Sts.Actions.*;
import com.Sts.Help.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.StsMath;

import javax.swing.*;
import java.awt.*;
import java.io.*;
import java.util.*;

public class Sts3dMenuBar extends JMenuBar
{
    JMenu menuFile = new JMenu();
    public JMenu menuToolbars = new JMenu();

        StsMenuItem menuFileNew = new StsMenuItem();
            JMenu menuNewProjects = new JMenu("Previous projects");
            StsMenuItem newProjectMenuItem = new StsMenuItem();
        StsMenuItem menuFileOpen = new StsMenuItem();
            JMenu menuDatabases = new JMenu("Previous databases");
            StsMenuItem otherDatabaseMenuItem = new StsMenuItem();

//        StsMenuItem menuWorkflow = new StsMenuItem();

//        StsMenuItem menuFileSave = new StsMenuItem();
        StsMenuItem menuFileSaveAs = new StsMenuItem();
        StsMenuItem menuFileExit = new StsMenuItem();

        JMenu menuArchiveProjects = new JMenu();
        StsMenuItem menuFileArchive = new StsMenuItem();
        StsMenuItem menuFileExplodeArchive = new StsMenuItem();

//    JMenu menuPreferences = new JMenu();
//        StsMenuItem menuPreferColor = new StsMenuItem();

//        StsMenuItem menuFileSetRoot = new StsMenuItem();
/*
        JMenu menuAsciiRead = new JMenu();
            StsMenuItem menuReadOWFaultPolygons = new StsMenuItem();
            StsMenuItem menuReadOWFaultCuts = new StsMenuItem();

            StsMenuItem menuReadWellPath = new StsMenuItem();
            StsMenuItem menuReadSurface = new StsMenuItem();
            StsMenuItem menuReadWellStratZone = new StsMenuItem();
            StsMenuItem menuReadLogCurve = new StsMenuItem();

        JMenu menuHDFRead = new JMenu();
            StsMenuItem menuReadHDFSurface = new StsMenuItem();
            StsMenuItem menuReadAndDecimateHDFSurface = new StsMenuItem();
            StsMenuItem menuReadHDFWellCatalog = new StsMenuItem();
            StsMenuItem menuReadHDFWell = new StsMenuItem();
		    StsMenuItem menuReadHDFPropertyVolume = new StsMenuItem();
        JMenu menuExport = new JMenu();
            StsMenuItem menuZoneProperties = new StsMenuItem();
            StsMenuItem menuZoneCategoricalProperties = new StsMenuItem();
            StsMenuItem menuExportToEclipse = new StsMenuItem();
            StsMenuItem menuExportOWFaultCuts = new StsMenuItem();
//            StsMenuItem menuExportToVIP = new StsMenuItem();

    JMenu menuBuild = new JMenu();
        StsMenuItem menuHorizonsWizard = new StsMenuItem();
        StsMenuItem menuZonesWizard = new StsMenuItem();
        JMenu menuFaults = new JMenu();
            StsMenuItem menuAddFaultOnCursor = new StsMenuItem();
            StsMenuItem menuAddFaultOnSection = new StsMenuItem();
            StsMenuItem menuAddSection = new StsMenuItem();
            StsMenuItem menuAddSectionEdgeOnSurface = new StsMenuItem();
            StsMenuItem menuAddEdgeOnSurface = new StsMenuItem();
            StsMenuItem menuEdgesFromImportedEdges = new StsMenuItem();
            StsMenuItem menuIntersectSectionEdges = new StsMenuItem();
            StsMenuItem menuAddAuxiliaryEdgeOnSurface = new StsMenuItem();
            StsMenuItem menuBuildFrame = new StsMenuItem();
            StsCheckboxMenuItem checkboxImportedEdges = new StsCheckboxMenuItem();

        JMenu menuInterfaces = new JMenu();
            StsMenuItem menuMarkersFromSurfaces = new StsMenuItem();
            StsMenuItem menuReorderMarkers = new StsMenuItem();
            StsMenuItem menuSurfacesFromMarkers = new StsMenuItem();
            StsMenuItem menuSurfacesFromIsopachs = new StsMenuItem();
            StsMenuItem menuReorderSurfaces = new StsMenuItem();
            StsMenuItem menuLinkSurfacesAndMarkers = new StsMenuItem();

        JMenu menuIntervals = new JMenu();
            StsMenuItem menuBuildZones = new StsMenuItem();
            StsMenuItem menuBuildWellZones = new StsMenuItem();

        JMenu menuBoundaries = new JMenu();
            StsMenuItem menuAddBoundary = new StsMenuItem();
            StsMenuItem menuAddRowColBoundary = new StsMenuItem();
            StsMenuItem menuBuildModelSurfaces = new StsMenuItem();

        JMenu menuModel = new JMenu();
            StsMenuItem menuBlockGrids = new StsMenuItem();
            StsMenuItem menuBuildModel = new StsMenuItem();

        JMenu menuWells = new JMenu();
            StsMenuItem menuAddPseudoOnCursor = new StsMenuItem();
            StsMenuItem menuAddPseudoOnSection = new StsMenuItem();

    JMenu menuEdit = new JMenu();
        StsMenuItem menuEditUndo = new StsMenuItem();

        JMenu menuDelete = new JMenu();
            StsMenuItem menuDeleteWell = new StsMenuItem();
            StsMenuItem menuDeleteSection = new StsMenuItem();
            StsMenuItem menuDeleteFaultPolygons = new StsMenuItem();

        JMenu menuEditColors = new JMenu();
            StsMenuItem menuToolsColor = new StsMenuItem();
            StsMenuItem menuEditSurfaceColors = new StsMenuItem();
            StsMenuItem menuEditZoneColors = new StsMenuItem();
            StsMenuItem menuEditWellColors = new StsMenuItem();
            StsMenuItem menuEditWellZoneColors = new StsMenuItem();
            StsMenuItem menuEditMarkerColors = new StsMenuItem();
            StsMenuItem menuEditLogCurveColors = new StsMenuItem();

        JMenu menuEditNames = new JMenu();
            StsMenuItem menuEditZoneSurfaceNames = new StsMenuItem();
            StsMenuItem menuEditZoneNames = new StsMenuItem();
            StsMenuItem menuEditWellZoneNames = new StsMenuItem();
*/
    JMenu menuView = new JMenu();
        // menuToolbars are added here after toolbars are generated
//            StsCheckboxMenuItem menuDefault = new StsCheckboxMenuItem();
//            StsCheckboxMenuItem menuViewSelect = new StsCheckboxMenuItem();
//            StsCheckboxMenuItem menuMouseControl = new StsCheckboxMenuItem();
//            StsCheckboxMenuItem menuCollaboration = new StsCheckboxMenuItem();
//            StsCheckboxMenuItem menuIntraFamily = new StsCheckboxMenuItem();
//            StsCheckboxMenuItem objectSelection = new StsCheckboxMenuItem();

        JMenu menuWindowControls = new JMenu();
            StsCheckboxMenuItem menuVisibilityDetail = new StsCheckboxMenuItem();
            StsCheckboxMenuItem menuViewProperties = new StsCheckboxMenuItem();
//            StsCheckboxMenuItem checkboxView3dCursors = new StsCheckboxMenuItem();

        JMenu menu3dGraphics = new JMenu();
            //StsCheckboxMenuItem checkboxPerspective = new StsCheckboxMenuItem();
            StsCheckboxMenuItem checkboxUseDisplayLists = new StsCheckboxMenuItem();
            StsCheckboxMenuItem checkboxViewShift = new StsCheckboxMenuItem();
            StsCheckboxMenuItem checkboxDrawConvex = new StsCheckboxMenuItem();

        JMenu menuDebug = new JMenu();
            StsMenuItem debugPicker = new StsMenuItem();
            StsMenuItem debugSurfacePoint = new StsMenuItem();
 //           StsMenuItem debugBuildProfiles = new StsMenuItem();
            StsCheckboxMenuItem checkboxConcavePolygonTest = new StsCheckboxMenuItem();
//            StsCheckboxMenuItem checkboxDebugGLDraw = new StsCheckboxMenuItem();
            StsCheckboxMenuItem checkboxDebugGridWeights = new StsCheckboxMenuItem();
            StsCheckboxMenuItem checkboxDebugSurfaceGapPoints = new StsCheckboxMenuItem();
            StsCheckboxMenuItem checkboxDebugBlockGapPoints = new StsCheckboxMenuItem();
            StsCheckboxMenuItem checkboxDebugNearestFaultIndices = new StsCheckboxMenuItem();
            StsCheckboxMenuItem checkboxDebugBlockZoneSidePoints = new StsCheckboxMenuItem();
            StsCheckboxMenuItem checkboxDebugRowLinksID = new StsCheckboxMenuItem();
            StsCheckboxMenuItem checkboxDebugColLinksID = new StsCheckboxMenuItem();
            StsCheckboxMenuItem checkboxDebugRowLinksSeq = new StsCheckboxMenuItem();
            StsCheckboxMenuItem checkboxDebugColLinksSeq = new StsCheckboxMenuItem();
            StsCheckboxMenuItem checkboxDebugDisplayRows = new StsCheckboxMenuItem();
            StsCheckboxMenuItem checkboxDebugDisplayCols = new StsCheckboxMenuItem();
            StsCheckboxMenuItem checkboxDebugDisplayBlockGrids = new StsCheckboxMenuItem();
            StsCheckboxMenuItem checkboxDebugDomainMap = new StsCheckboxMenuItem();
            StsCheckboxMenuItem checkboxDebugGridDomains = new StsCheckboxMenuItem();
            StsCheckboxMenuItem checkboxDebugDisplayZoneSides = new StsCheckboxMenuItem();
            StsCheckboxMenuItem checkboxDebugFrame = new StsCheckboxMenuItem();
            StsCheckboxMenuItem checkboxDebugModel = new StsCheckboxMenuItem();
            StsCheckboxMenuItem checkboxDebugExceptions = new StsCheckboxMenuItem();
/*
        JMenu menuListings = new JMenu();
            //StsMenuItem menuListSurfaces = new StsMenuItem();
            StsMenuItem menuListRawSurfaces = new StsMenuItem();
            StsMenuItem menuListZoneSurfaces = new StsMenuItem();
            StsMenuItem menuListZones = new StsMenuItem();
            StsMenuItem menuListWells = new StsMenuItem();
            StsMenuItem menuListWellZones = new StsMenuItem();
            StsMenuItem menuListWellZoneSets = new StsMenuItem();
            StsMenuItem menuListWellMarkers = new StsMenuItem();
            StsMenuItem menuListMarkers = new StsMenuItem();
            StsMenuItem menuListLogs = new StsMenuItem();
            StsMenuItem menuPrintLogs = new StsMenuItem();
            StsMenuItem menuListLogSets = new StsMenuItem();
*/
        StsMenuItem menuReportSurface = new StsMenuItem();

//        StsMenuItem menuDisplayLogFile = new StsMenuItem();

//        JMenu menuProperties = new JMenu();
//            StsMenuItem menuPropertyOnSurface= new StsMenuItem();
//            StsMenuItem menuRemoveSurfaceProperty= new StsMenuItem();

    JMenu menuHelp = new JMenu();
        StsMenuItem menuHelpAbout = new StsMenuItem();
        StsMenuItem menuUserManual = new StsMenuItem();
        StsMenuItem menuStepsManual = new StsMenuItem();
        StsMenuItem menuUsageInfo = new StsMenuItem();

	private StsModel model;
	private StsWin3d win3d;
//	private StsGLPanel3d glPanel3d;
	private StsActionManager windowActionManager;

	public Sts3dMenuBar(StsModel model, StsWin3d win3d)
	{
        this.model = model;
		this.win3d = win3d;
//    	this.glPanel3d = glPanel3d;
		this.windowActionManager = win3d.getActionManager();
       	try
        {
        	jbInit();
       	}
        catch (Exception e)
        {
            e.printStackTrace();
        }
	}

//Component initialization

    private void jbInit() throws Exception
    {
        menuFile.setText("File");
        add(menuFile);
        buildMenuFile();

//        menuPreferences.setLabel("Preferences");
//        add(menuPreferences);
//        buildMenuPreferences();
/*
           	menuFile.add(menuHDFRead);
			menuHDFRead.setLabel("Read Sts Data ");
	    	    menuHDFRead.add(menuReadHDFSurface);
			    menuReadHDFSurface.setMenuAction("Surfaces...",
                        StsImportHDFSurfaces.class, actionManager);

	    	    menuHDFRead.add(menuReadAndDecimateHDFSurface);
			    menuReadAndDecimateHDFSurface.setMenuAction("Decimated Surfaces...",
                        StsImportAndDecimateHDFSurfaces.class, actionManager);

	    	    menuHDFRead.add(menuReadHDFWellCatalog);
			    menuReadHDFWellCatalog.setMenuAction("Wells from Catalog...",
                        StsImportHDFWellsFromCatalog.class, actionManager);

	    	    menuHDFRead.add(menuReadHDFWell);
			    menuReadHDFWell.setMenuAction("Wells from Files...",
                        StsImportHDFWells.class, actionManager);

                menuHDFRead.addSeparator();

	    	    menuHDFRead.add(menuReadHDFPropertyVolume);
			    menuReadHDFPropertyVolume.setMenuAction("3-D Model Properties...",
                        StsImportHDFPropertyVolume.class, actionManager);

            menuFile.add(menuExport);
            menuExport.setLabel("Export ");
           		menuExport.add(menuZoneProperties);
                menuZoneProperties.setMenuAction(
                        "Continuous Log Properties for ResMod...",
                        StsExportZoneLogProperties.class, actionManager);

           		menuExport.add(menuZoneCategoricalProperties);
                menuZoneCategoricalProperties.setMenuAction(
                        "Categorical Log Properties for ResMod...",
                        StsExportCategoricalFacies.class, actionManager);

                menuExport.addSeparator();

                menuExport.add(menuExportOWFaultCuts);
                menuExportOWFaultCuts.setMenuAction("Fault Cuts...",
                        StsExportFaultCuts.class, actionManager);

                menuExport.addSeparator();

                // !! replace this with a real action
                menuExport.add(menuExportToEclipse);
        		menuExportToEclipse.setMenuActionListener(
                        "3-D Grid to Eclipse...",
                        glPanel, "statusBarMessage",
                        new String("Export to Eclipse not yet implemented!"));
*/
/*
                // !! replace this with a real action
                menuExport.add(menuExportToVIP);
        		menuExportToVIP.setMenuActionListener(
                        "3-D Grid to VIP...",
                        glPanel3d, "statusBarMessage",
                        new String("Export to VIP not yet implemented!"));
*/
//            menuFile.addSeparator();
/*
           	menuFile.add(menuAsciiRead);
			menuAsciiRead.setLabel("Import Ascii ");
                menuAsciiRead.add(menuReadOWFaultPolygons);
                menuReadOWFaultPolygons.setMenuAction("OpenWorks Fault Polygons...",
                        StsImportOWFaultPolygons.class, actionManager);

                menuAsciiRead.add(menuReadOWFaultCuts);
                menuReadOWFaultCuts.setMenuAction("OpenWorks Fault Cuts...",
                        StsImportFaultCuts.class, actionManager);

	    	    menuAsciiRead.add(menuReadSurface);
		        menuReadSurface.setMenuAction("Map Surfaces...",
                        StsImportAsciiSurfaces.class, actionManager);

                menuAsciiRead.addSeparator();

	    	    menuAsciiRead.add(menuReadWellPath);
			    menuReadWellPath.setMenuAction("Well Deviations...",
                        StsImportAsciiWells.class, actionManager);


	    	    menuAsciiRead.add(menuReadWellStratZone);
		        menuReadWellStratZone.setMenuAction("Well Strat Zones...",
                        StsImportAsciiWellStratZones.class, actionManager);

	    	    menuAsciiRead.add(menuReadLogCurve);
		        menuReadLogCurve.setMenuAction("Multiple Log Curves...",
                        StsImportAsciiLogCurves.class, actionManager);

			menuFile.add(menuFileSetRoot);
        	menuFileSetRoot.setMenuAction("Set Root Directory...",
                    StsSetRootDirectory.class, actionManager);

            menuFile.addSeparator();
*/
//			menuFile.add(menuFileExit);
 //       	menuFileExit.setMenuActionListener("Exit", win3d, "appClose", null);
/*
		add(menuBuild);
		menuBuild.setLabel("Build");
		    menuBuild.add(menuInterfaces);
                menuInterfaces.setLabel("Horizons ");
	    	    menuInterfaces.add(menuHorizonsWizard);
                menuHorizonsWizard.setMenuAction("Horizon Wizard...",
                            StsHorizonsWizard.class, actionManager);

                menuInterfaces.addSeparator();

                menuInterfaces.add(menuMarkersFromSurfaces);
                menuMarkersFromSurfaces.setMenuAction("Markers From Surfaces...",
                        StsBuildMarkersFromSurfaces.class, actionManager);

                menuInterfaces.add(menuReorderMarkers);
                menuReorderMarkers.setMenuAction("Reorder Markers...",
                        StsReorderMarkers.class, actionManager);

                menuInterfaces.addSeparator();

                menuInterfaces.add(menuSurfacesFromMarkers);
                menuSurfacesFromMarkers.setMenuAction("Surfaces From Markers...",
                        StsBuildSurfacesFromMarkers.class, actionManager);

                menuInterfaces.add(menuSurfacesFromIsopachs);
                menuSurfacesFromIsopachs.setMenuAction("Surfaces From Isopachs...",
                        StsBuildSurfacesFromIsopachs.class, actionManager);

     			menuInterfaces.add(menuReorderSurfaces);
            	menuReorderSurfaces.setMenuAction("Reorder Surfaces...",
                        StsReorderSurfaces.class, actionManager);

                menuInterfaces.addSeparator();

                menuInterfaces.add(menuLinkSurfacesAndMarkers);
                menuLinkSurfacesAndMarkers.setMenuAction("Correlate Surfaces And Markers...",
                        StsLinkSurfaceWithMarker.class, actionManager);

		    menuBuild.add(menuIntervals);
            menuIntervals.setLabel("Zones ");
    		    menuIntervals.add(menuZonesWizard);
                menuZonesWizard.setMenuAction("Zone Wizard...",
                            StsModelWizard.class, actionManager);

                menuIntervals.addSeparator();

			    menuIntervals.add(menuBuildZones);
            	menuBuildZones.setMenuAction("Select Horizons for Zones...",
                        StsBuildZonesFromSurfaces.class, actionManager);

                menuIntervals.add(menuBuildModelSurfaces);
                menuBuildModelSurfaces.setMenuAction("Build Horizon Grids",
                        StsBuildModelSurfaces.class, actionManager);

		    menuBuild.add(menuBoundaries);
            menuBoundaries.setLabel("Boundaries ");
                menuBoundaries.add(menuAddBoundary);
                menuAddBoundary.setMenuAction("Pick Arbitrary Boundary...",
                        StsPolygonBoundary.class, actionManager);

                menuBoundaries.add(menuAddRowColBoundary);
                menuAddRowColBoundary.setMenuAction("Pick Row-Col Boundary...",
                        StsRectangularBoundary.class, actionManager);


		    menuBuild.add(menuFaults);
            menuFaults.setLabel("Faults");
                if (Main.NON_VERTICAL_FAULTS)
                {
                    menuFaults.add(menuAddFaultOnCursor);
                    menuAddFaultOnCursor.setMenuAction("Fault Line on Cursor...",
                            StsWellOnCursor.class, actionManager);

                    menuFaults.add(menuAddFaultOnSection);
                    menuAddFaultOnSection.setMenuAction("Fault Line on Section...",
                            StsWellOnSection.class, actionManager);

                    menuFaults.add(menuAddSection);
                    menuAddSection.setMenuAction("Fault Section Between Fault Lines...",
                            StsAddSection.class, actionManager);

                    menuFaults.add(menuAddSectionEdgeOnSurface);
                    menuAddSectionEdgeOnSurface.setMenuAction("Surface Cut Between Fault Lines...",
                            StsSectionEdgeOnSurface.class, actionManager,
                            new Object[]{new Integer(StsSection.FAULT)});

                    menuFaults.addSeparator();
                }

                menuFaults.add(menuAddEdgeOnSurface);
        	    args = new Object[] { new Integer(StsParameters.FAULT) };
                menuAddEdgeOnSurface.setMenuAction("Surface Fault Cuts...",
                        StsEdgeOnSurface.class, actionManager, args);

                menuFaults.addSeparator();

                menuFaults.add(menuEdgesFromImportedEdges);
                menuEdgesFromImportedEdges.setMenuAction("Convert imported fault edges to model edges",
                        StsEdgesFromImportedEdges.class, actionManager);

                menuFaults.add(menuIntersectSectionEdges);
                menuIntersectSectionEdges.setMenuAction("Automatically intersect/split sectionEdges",
                        StsIntersectSectionEdges.class, actionManager);

                menuFaults.add(menuAddAuxiliaryEdgeOnSurface);
        	    args = new Object[] { new Integer(StsParameters.AUXILIARY) };
                menuAddAuxiliaryEdgeOnSurface.setMenuAction("Surface Auxiliary Cuts...",
                        StsEdgeOnSurface.class, actionManager, args);


                menuFaults.add(menuBuildFrame);
                menuBuildFrame.setMenuAction("Fault Sections From Surface Cuts",
                        StsBuildSections.class, actionManager);

                menuFaults.add(checkboxImportedEdges);
                checkboxImportedEdges.setMenuItemListener("Toggle imported fault edges",
                        model, "toggleImportedFaults");

                checkboxImportedEdges.setState(true);


		    menuBuild.add(menuModel);
            menuModel.setLabel("Model ");
                menuModel.add(menuBlockGrids);
                menuBlockGrids.setMenuAction("Test Fault Gapping",
                        StsGapBlockGrids.class, actionManager);

                menuModel.addSeparator();

                menuModel.add(menuBuildModel);
                menuBuildModel.setMenuAction("Build 3-D Zone Grids",
                        StsBuildModel.class, actionManager);

            if (Main.NON_VERTICAL_FAULTS)
            {
                menuBuild.addSeparator();

                menuBuild.add(menuWells);
                menuWells.setLabel("Wells ");
                    menuWells.add(menuAddPseudoOnCursor);
                    menuAddPseudoOnCursor.setMenuAction("Pseudo on Cursor...",
                            StsWellOnCursor.class, actionManager);

                    menuWells.add(menuAddPseudoOnSection);
                    menuAddPseudoOnSection.setMenuAction("Pseudo on Section...",
                            StsWellOnSection.class, actionManager);
            }

		add(menuEdit);
		menuEdit.setLabel("Edit");

			menuEdit.add(menuEditUndo);
			menuEditUndo.setMenuActionListener("Undo", StsWindowActionManager.class, "undoAction", null);

    		menuEdit.add(menuDelete);
            menuDelete.setLabel("Delete ");
                menuDelete.add(menuDeleteWell);
                menuDeleteWell.setMenuAction("Well, Fault Line, or Pseudo-Well...",
                        StsDeleteWell.class, actionManager);

                menuDelete.add(menuDeleteSection);
                menuDeleteSection.setMenuAction("Section...",
                        StsDeleteSection.class, actionManager);

                menuDelete.add(menuDeleteFaultPolygons);
                menuDeleteFaultPolygons.setMenuAction("Fault Polygons...",
                        StsDeleteFaultPolygons.class, actionManager);

            menuEdit.addSeparator();

            menuEdit.add(menuEditColors);
            menuEditColors.setLabel("Change Colors ");
    			menuEditColors.add(menuToolsColor);
 	    		menuToolsColor.setMenuActionListener("Current Color...",
                        model, "displaySpectrum", new String("Basic"));

                menuEditColors.addSeparator();

                menuEditColors.add(menuEditSurfaceColors);
        		menuEditSurfaceColors.setMenuAction("Surfaces & Horizons...",
                        StsEditSurfaceColors.class, actionManager);

                menuEditColors.add(menuEditZoneColors);
        		menuEditZoneColors.setMenuAction("Zones...",
                        StsEditZoneColors.class, actionManager);

                menuEditColors.addSeparator();

                menuEditColors.add(menuEditWellColors);
        		menuEditWellColors.setMenuAction("Wells...",
                        StsEditWellColors.class, actionManager);

                menuEditColors.add(menuEditWellZoneColors);
        		menuEditWellZoneColors.setMenuAction("Well Zones...",
                        StsEditWellZoneColors.class, actionManager);

                menuEditColors.add(menuEditMarkerColors);
        		menuEditMarkerColors.setMenuAction("Well Markers...",
                        StsEditMarkerColors.class, actionManager);

                menuEditColors.add(menuEditLogCurveColors);
        		menuEditLogCurveColors.setMenuAction("Log Curves...",
                        StsEditLogCurveColors.class, actionManager);

            menuEdit.add(menuEditNames);
            menuEditNames.setLabel("Rename ");
                menuEditNames.add(menuEditZoneSurfaceNames);
        		menuEditZoneSurfaceNames.setMenuAction("Horizons...",
                        StsRenameZoneSurfaces.class, actionManager);

                menuEditNames.add(menuEditZoneNames);
        		menuEditZoneNames.setMenuAction("Zones...",
                        StsRenameZones.class, actionManager);

                menuEditNames.addSeparator();

                menuEditNames.add(menuEditWellZoneNames);
        		menuEditWellZoneNames.setMenuAction("Well Zones...",
                        StsRenameWellZones.class, actionManager);
*/
 		add(menuView);
		menuView.setText("View");
        menuView.add(menuToolbars);
            menuToolbars.setText("Toolbars");
        /*
            win3d.addToolbarMenuItems(menuToolbars);
            menuToolbars.add(win3d.menubarDefault);
			win3d.menubarDefault.setState(true);
			win3d.menuDefault.setState(true);
            win3d.menubarDefault.setMenuItemListener("Default", win3d, "defaultToolbar");
            win3d.menuDefault.setMenuItemListener("Default", win3d, "defaultToolbar");

            menuToolbars.add(win3d.menubarViewSelect);
			win3d.menubarViewSelect.setState(true);
			win3d.menuViewSelect.setState(true);
            win3d.menubarViewSelect.setMenuItemListener("View Select", win3d, "viewSelectToolbar");
            win3d.menuViewSelect.setMenuItemListener("View Select", win3d, "viewSelectToolbar");

            menuToolbars.add(win3d.menubarMouseControl);
			win3d.menubarMouseControl.setState(true);
			win3d.menuMouseControl.setState(true);
            win3d.menubarMouseControl.setMenuItemListener("Mouse Control", win3d, "mouseActionToolbar");
            win3d.menuMouseControl.setMenuItemListener("Mouse Control", win3d, "mouseActionToolbar");


            menuToolbars.add(win3d.menubarTime);
            win3d.menubarTime.setState(false);
            win3d.menuTime.setState(false);
            win3d.menubarTime.setMenuItemListener("Time Control", win3d, "timeActionToolbar");
            win3d.menuTime.setMenuItemListener("Time Control", win3d, "timeActionToolbar");

            menuToolbars.add(win3d.menubarIntraFamily);
			win3d.menubarIntraFamily.setState(true);
			win3d.menuIntraFamily.setState(true);
            win3d.menubarIntraFamily.setMenuItemListener("IntraFamily", win3d, "intraFamilyToolbar");
            win3d.menuIntraFamily.setMenuItemListener("IntraFamily", win3d, "intraFamilyToolbar");

            menuToolbars.add(win3d.menubarMedia);
			win3d.menubarMedia.setState(false);
			win3d.menuMedia.setState(false);
            win3d.menubarMedia.setMenuItemListener("Media", win3d, "mediaToolbar");
            win3d.menuMedia.setMenuItemListener("Media", win3d, "mediaToolbar");

            menuToolbars.add(win3d.menubarObjectSelection);
			win3d.menubarObjectSelection.setState(true);
			win3d.objectSelection.setState(true);
            win3d.menubarObjectSelection.setMenuItemListener("Object Selection", win3d, "comboBoxToolbar");
            win3d.objectSelection.setMenuItemListener("Object Selection", win3d, "comboBoxToolbar");

            menuToolbars.add(win3d.menubarCollaboration);
            win3d.menubarCollaboration.setState(true);
            win3d.menuCollaboration.setState(true);
            win3d.menubarCollaboration.setMenuItemListener("Collaboration", win3d, "collaborationToolbar");
            win3d.menuCollaboration.setMenuItemListener("Collaboration", win3d, "collaborationToolbar");
        */
        menuView.add(menuWindowControls);
            menuWindowControls.setText("Window Controls ");
            menuWindowControls.add(menuVisibilityDetail);
			menuVisibilityDetail.setState(true);
            menuVisibilityDetail.setMenuItemListener("Show Tree View", win3d, "toggleVisibilityPane");

//            menuWindowControls.add(menuViewProperties);
//            menuViewProperties.setMenuItemListener("Show 3-D Controls", win3d, "toggleView3dPane");
//            menuViewProperties.setState(true);
//
//        	    menuWindowControls.add(checkboxView3dCursors);
//                checkboxView3dCursors.setMenuItemListener("Sbow 3-D Cursor",
//                        win3d, "toggleView3dCursors");

        menuView.add(menu3dGraphics);
            menu3dGraphics.setText("Display Options ");

            // Removed from menubar and added to project properties, don't need in both places.
            //menu3dGraphics.add(this.checkboxPerspective);
            //checkboxPerspective.setState(win3d.isPerspective);
            //checkboxPerspective.setMenuItemListener("Perspective Projection", win3d, "toggleProjection");

            menu3dGraphics.add(checkboxUseDisplayLists);
			checkboxUseDisplayLists.setState(true);
			checkboxUseDisplayLists.setToggleMenuItemListener("Use Display Lists", model, "toggleUseDisplayLists");

     	    menu3dGraphics.add(checkboxViewShift);
		    checkboxViewShift.setState(true);
            checkboxViewShift.setMenuItemListener("Use View Offsets", win3d, "toggleViewShift");

      	    menu3dGraphics.add(checkboxDrawConvex);
			checkboxDrawConvex.setState(true);
            checkboxDrawConvex.setMenuItemListener("Draw Convex Only", StsPolygon.class, "toggleDrawConvex");

        menuView.addSeparator();

		menuView.add(menuDebug);
			menuDebug.setText("Debug ");
                menuDebug.add(debugPicker);
                debugPicker.setMenuAction("Debug picking", StsDebugPicker.class, windowActionManager);

                menuDebug.add(debugSurfacePoint);
//                debugSurfacePoint.setMenuAction("Surface Point mainDebug...", StsSurfacePointDebug.class, actionManager);
                debugSurfacePoint.setMenuAction("Surface Point mainDebug...", StsSurfacePointDebug.class, windowActionManager);
//                menuDebug.add(debugBuildProfiles);
//                debugSurfacePoint.setMenuAction("Surface Point mainDebug...", StsSurfacePointDebug.class, actionManager);
//                debugSurfacePoint.setMenuAction("Debug build profiles", StsDebugBuildProfiles.class, actionManager);

                menuDebug.add(checkboxConcavePolygonTest);
				checkboxConcavePolygonTest.setState(false);
                checkboxConcavePolygonTest.setMenuItemListener("Concave Polygon Test", StsCursor3d.class, "toggleConcavePolygonTest");

//				menuDebug.add(checkboxDebugGLDraw);
//				checkboxDebugGLDraw.setState(false);
//                checkboxDebugGLDraw.setMenuItemListener("GLDraw", StsGLPanel.class, "toggleGLDrawDebug");

				menuDebug.add(checkboxDebugGridWeights);
				checkboxDebugGridWeights.setState(false);
                checkboxDebugGridWeights.setToggleMenuItemListener("Grid Weights", model, "togglePropertyDisplayAll");

				menuDebug.add(checkboxDebugBlockGapPoints);
				checkboxDebugBlockGapPoints.setState(false);
				checkboxDebugBlockGapPoints.setToggleMenuItemListener("Block Gap Points", model, "togglePropertyDisplayAll");

                menuDebug.add(checkboxDebugSurfaceGapPoints);
				checkboxDebugSurfaceGapPoints.setState(false);
				checkboxDebugSurfaceGapPoints.setToggleMenuItemListener("Surface Gap Points", model, "togglePropertyDisplayAll");
          /*
				menuDebug.add(checkboxDebugNearestFaultIndices);
                checkboxDebugNearestFaultIndices.setMenuItemListener("Nearest Fault Indices",
                        model, "toggleDebugNearestFaultIndices");
                checkboxDebugNearestFaultIndices.setState(false);
          */
				menuDebug.add(checkboxDebugBlockZoneSidePoints);
				checkboxDebugBlockZoneSidePoints.setState(false);
                checkboxDebugBlockZoneSidePoints.setToggleMenuItemListener("BlockZoneSide Points", model, "togglePropertyDisplayAll");

				menuDebug.add(checkboxDebugRowLinksID);
				checkboxDebugRowLinksID.setState(false);
                checkboxDebugRowLinksID.setToggleMenuItemListener("Row Links ID", model, "togglePropertyDisplayAll");

				menuDebug.add(checkboxDebugColLinksID);
				checkboxDebugColLinksID.setState(false);
                checkboxDebugColLinksID.setToggleMenuItemListener("Column Links ID", model, "togglePropertyDisplayAll");

				menuDebug.add(checkboxDebugRowLinksSeq);
				checkboxDebugRowLinksSeq.setState(false);
                checkboxDebugRowLinksSeq.setToggleMenuItemListener("Row Links Seq", model, "togglePropertyDisplayAll");

				menuDebug.add(checkboxDebugColLinksSeq);
				checkboxDebugColLinksSeq.setState(false);
                checkboxDebugColLinksSeq.setToggleMenuItemListener("Column Links Seq", model, "togglePropertyDisplayAll");

				menuDebug.add(checkboxDebugDisplayRows);
				checkboxDebugDisplayRows.setState(false);
                checkboxDebugDisplayRows.setToggleMenuItemListener("Display Rows", model, "togglePropertyDisplayAll");

				menuDebug.add(checkboxDebugDisplayCols);
				checkboxDebugDisplayCols.setState(false);
                checkboxDebugDisplayCols.setToggleMenuItemListener("Display Columns", model, "togglePropertyDisplayAll");

				menuDebug.add(checkboxDebugDisplayBlockGrids);
				checkboxDebugDisplayBlockGrids.setState(false);
                checkboxDebugDisplayBlockGrids.setToggleMenuItemListener("Display Gap Grids", model, "togglePropertyDisplayAll");

				menuDebug.add(checkboxDebugGridDomains);
				checkboxDebugGridDomains.setState(false);
                checkboxDebugGridDomains.setToggleMenuItemListener("Display Cell Types", model, "togglePropertyDisplayAll");

				menuDebug.add(checkboxDebugDisplayZoneSides);
				checkboxDebugDisplayZoneSides.setState(false);
                checkboxDebugDisplayZoneSides.setToggleMenuItemListener("Display Zone Sides", model, "togglePropertyDisplayAll");

				menuDebug.add(checkboxDebugFrame);
				checkboxDebugFrame.setState(false);
                checkboxDebugFrame.setToggleMenuItemListener("Frame", model, "toggleProperty");

				menuDebug.add(checkboxDebugModel);
				checkboxDebugModel.setState(false);
                checkboxDebugModel.setToggleMenuItemListener("Model", model, "toggleProperty");

				menuDebug.add(checkboxDebugExceptions);
				checkboxDebugExceptions.setState(false);
                checkboxDebugExceptions.setToggleMenuItemListener("Exceptions", model, "toggleProperty");

            menuView.add(menuReportSurface);
            menuReportSurface.setMenuAction("Surface Reporting...", StsSurfaceReporting.class, windowActionManager);

            menuView.addSeparator();

		add(menuHelp);
        menuHelp.setText("Help");
 			menuHelp.add(menuHelpAbout);
       		menuHelpAbout.setMenuActionListener("About...", win3d, "helpAboutDialog", null);

 			menuHelp.add(menuUserManual);
            menuHelp.add(menuStepsManual);

            menuUserManual.setLabel("General Operations Manual...");
            HelpManager.setSourceHelpActionListener(menuUserManual, HelpManager.GENERAL, model.win3d);  // JavaHelp

            menuStepsManual.setLabel("Workflow Steps Manual...");
            HelpManager.setSourceHelpActionListener(menuStepsManual, HelpManager.STEPS, model.win3d);  // JavaHelp

            menuHelp.add(menuUsageInfo);
       		menuUsageInfo.setMenuActionListener("Usage...", win3d, "usageDialog", null);
	}
/*
    public void loadManual(StsMenuItem m, String name)
    {
        HelpManager.getRemoteHelp(m, name + ".hs", name + "Help.jar", model.win3d);
    }
*/

    private void buildMenuFile()
    {
        menuFile.removeAll();

        String[] projectDirectories = model.getProject().getProjectDirectories();
        int nProjects = projectDirectories.length;
        if(nProjects == 0)
        {
            menuFileNew.setMenuAction("New", StsNewModel.class, windowActionManager);
            menuFile.add(menuFileNew);
        }
        else
        {
            menuNewProjects = new JMenu("New");
            JMenuItem menuLabel = new JMenuItem("Select Project Directory from List Below");
            menuLabel.setFont(new Font("Dialog", Font.BOLD, 12));
            menuNewProjects.add(menuLabel);
            menuNewProjects.addSeparator();

            // Sort into last modified
            dbTimeSort[] dbs = new dbTimeSort[nProjects];
            for(int n=0; n< nProjects; n++)
            {
                File tFile = new File(projectDirectories[n]);
                dbs[n] = new dbTimeSort(projectDirectories[n], tFile.lastModified());
            }
            Arrays.sort(dbs);

            StsMenuItem projectMenuItem;
            nProjects = Math.min(nProjects, 20);
            for(int n = 0; n < nProjects; n++)
            {
                projectMenuItem = new StsMenuItem();
                projectMenuItem.setMenuAction(dbs[n].dbName, StsNewModel.class, windowActionManager, new Object[] { dbs[n].dbName } );
                menuNewProjects.add(projectMenuItem);
            }
            menuNewProjects.addSeparator();

            newProjectMenuItem.setMenuAction("New...", StsNewModel.class, windowActionManager);
            menuNewProjects.add(newProjectMenuItem);
            menuFile.add(menuNewProjects);
        }

        String[] projectDatabases = StsProject.getDatabases();
        int nDatabases = projectDatabases.length;
        if(nDatabases == 0)
        {
            menuFileOpen.setMenuAction("Open", StsOpenModel.class, windowActionManager);
            menuFile.add(menuFileOpen);
        }
        else
        {
            menuDatabases = new JMenu("Open");
            buildMenuDatabases();
            menuFile.add(menuDatabases);
        }
        menuArchiveProjects = new JMenu("Archive");

        menuFileArchive.setMenuActionListener("Create...", model, "createArchive", null);
        menuFileExplodeArchive.setMenuActionListener("Explode...", model, "explodeArchive", windowActionManager);
        menuArchiveProjects.add(menuFileArchive);
        menuArchiveProjects.add(menuFileExplodeArchive);
        menuFile.add(menuArchiveProjects);

        menuFile.addSeparator();

        Object[] args;

//        args = new Object[] { new Integer(StsSaveModel.SAVE) };
//        menuFileSave.setMenuAction("Save", StsSaveModel.class, actionManager, args);
//        menuFile.add(menuFileSave);

        args = new Object[] { new Integer(StsSaveModel.SAVE_AS) };
        menuFileSaveAs.setMenuAction("Save As...", StsSaveModel.class, windowActionManager, args);
        menuFile.add(menuFileSaveAs);

        menuFile.addSeparator();

        menuFile.add(menuFileExit);
        menuFileExit.setMenuActionListener("Exit", model, "appClose", null);
    }

   public void buildMenuDatabases()
    {
        String[] projectDatabases = StsProject.getDatabases();
        int nDatabases = projectDatabases.length;

        // Sort into last modified

        dbTimeSort[] dbs = new dbTimeSort[nDatabases];
        int nExistingDatabases = 0;
        for(int n=0; n< nDatabases; n++)
        {
            File tFile = new File(new String(projectDatabases[n]));
            if(!tFile.exists()) continue;
            dbs[nExistingDatabases++] = new dbTimeSort(projectDatabases[n], tFile.lastModified());
        }
        dbs = (dbTimeSort[])StsMath.trimArray(dbs, nExistingDatabases);
        nDatabases = nExistingDatabases;
        Arrays.sort(dbs);

        if(menuDatabases.getItemCount() > 0)
            menuDatabases.removeAll();

        JMenuItem menuLabel = new JMenuItem("Select Database from List Below");
        menuLabel.setFont(new Font("Dialog", Font.BOLD, 12));
        menuDatabases.add(menuLabel);

        menuDatabases.addSeparator();

        StsMenuItem databaseMenuItem;
        nDatabases = Math.min(nDatabases, 20);
        for(int n = 0; n < nDatabases; n++)
        {
            databaseMenuItem = new StsMenuItem();
            databaseMenuItem.setMenuAction(dbs[n].dbName, StsOpenModel.class, windowActionManager, new Object[] { dbs[n].dbName } );
            menuDatabases.add(databaseMenuItem);
        }

        menuDatabases.addSeparator();

        otherDatabaseMenuItem.setMenuAction("Other...", StsOpenModel.class, windowActionManager);
        menuDatabases.add(otherDatabaseMenuItem);
    }

    /** sorts so that most recent object is first */
    class dbTimeSort implements Comparable
    {
        private String  dbName;
        private long lastMod;

        dbTimeSort(String dbName, long lastMod)
        {
            this.dbName = dbName;
            this.lastMod = lastMod;
        }

        public int compareTo(Object other)
        {
            if(lastMod < ((dbTimeSort)other).lastMod) return 1;
            return -1;
        }
    }

    public void addNewProject(String projectDirname)
    {
        model.getProject().addProjectDirectory(projectDirname);
    }

    /** enable/disable menu items */
    public void setEnabledMenuItems(boolean enable)
    {
//		menuFileSave.setEnabled(enable);
		menuFileSaveAs.setEnabled(enable);
 		menuView.setEnabled(enable);
    }

	public void setCheckboxUseDisplayLists(boolean enable)
    {
    	checkboxUseDisplayLists.setState(enable);
    }

    public JPopupMenu getToolbarPopupMenu()
    {
        JPopupMenu popupMenu = new JPopupMenu("Toolbars");
        StsCheckboxMenuItem[] popupToolbarItems = win3d.menubarToolbarItems;
        for(int n = 0; n < popupToolbarItems.length; n++)
            popupMenu.add(popupToolbarItems[n]);
        /*
        popupMenu.add(win3d.menuDefault);
        popupMenu.add(win3d.menuViewSelect);
        popupMenu.add(win3d.menuMouseControl);
        popupMenu.add(win3d.menuCollaboration);
        popupMenu.add(win3d.menuIntraFamily);
        popupMenu.add(win3d.menuMedia);
        popupMenu.add(win3d.menuTime);
        popupMenu.add(win3d.objectSelection);
        */
        return popupMenu;
    }
/*
    public void validateToolbars()
    {
        win3d.menuDefault.setState(win3d.defaultToolbar.isValid());
        win3d.menuViewSelect.setState(win3d.viewSelectToolbar.isValid());
        win3d.menuMouseControl.setState(win3d.mouseActionToolbar.isValid());
        win3d.menuIntraFamily.setState(win3d.intraFamilyToolbar.isValid());
        win3d.menuMedia.setState(win3d.mediaToolbar.isValid());
//        win3d.menuTime.setState(win3d.timeToolbar.isValid());
        win3d.objectSelection.setState(win3d.comboBoxToolbar.isValid());
 //       win3d.menubarCollaboration.setState(win3d.collaborationToolbar.isValid());
    }
*/
}