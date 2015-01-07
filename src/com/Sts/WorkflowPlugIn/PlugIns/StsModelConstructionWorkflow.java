package com.Sts.WorkflowPlugIn.PlugIns;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version 1.1
 */

import com.Sts.Workflow.*;
import com.Sts.WorkflowPlugIn.*;
import com.Sts.UI.ObjectPanel.StsObjectTree;
import com.Sts.UI.ObjectPanel.StsTreeNode;
import com.Sts.MVC.StsModel;
import com.Sts.DBTypes.StsVspClass;
import com.Sts.DBTypes.StsAncillaryDataClass;
import com.Sts.Types.StsSubType;
import com.Sts.Utilities.StsException;

public class StsModelConstructionWorkflow extends StsWorkflowPlugIn
{
    public StsModelConstructionWorkflow()
    {
        name = "StsModelConstructionWorkflow";
        workflowName = "Model Construction";
/*
        plugInClasses = new String[]
        {
        };
*/
        description = new String("Builds faulted model grid from horizons, faults and wells.");

        createComboBoxDescriptors();
    }

    public void createComboBoxDescriptors()
    {
        StsWorkflowPlugIn.ComboBoxDescriptor comboBoxDescriptor;

        String[] volumeClasses = new String[] {"com.Sts.DBTypes.StsSeismicVolume", "com.Sts.DBTypes.StsMathVirtualVolume", "com.Sts.DBTypes.StsBlendedVirtualVolume",
                "com.Sts.DBTypes.StsFilterVirtualVolume", "com.Sts.DBTypes.StsRGBAVirtualVolume",
                "com.Sts.DBTypes.StsCrossplotVirtualVolume", "com.Sts.DBTypes.StsSensorVirtualVolume" };
        comboBoxDescriptor = constructComboBoxToolbarDescriptor("com.Sts.DBTypes.StsSeismicVolume", volumeClasses, "seismic", "noseismic");
        comboBoxToolbarDescriptors.add(comboBoxDescriptor);
        comboBoxDescriptor = constructComboBoxDescriptor("com.Sts.DBTypes.StsCrossplot", "crossplot", "noCrossplot");
        comboBoxToolbarDescriptors.add(comboBoxDescriptor);
    }

    public void createWorkflowNodes(StsTreeModel treeModel, StsWorkflowTreeNode workflowRoot)
    {
        treeModel.addMenuNode(newProject);
        treeModel.addMenuNode(openProject);

        createGroupAndNodes(workflowRoot, treeModel, PROCESS_SEISMIC, new byte[] {P_POSTSTACK3D, P_POSTSTACK2D, P_VSP});

        createGroupAndNodes(workflowRoot, treeModel, LOAD_DATA, new byte[] {L_POSTSTACK3D, L_POSTSTACK2D, L_WELLSANDLOGS, L_VSP, L_SURFACE, L_FAULT_STICKS, L_COLORPALETTE, L_CULTURE, L_ANCILLARY});
        StsWorkflowTreeNode loadEclipseFile = getGroupNode(LOAD_DATA).addChild("com.Sts.Actions.Wizards.EclipseLoad.StsEclipseLoadWizard", "Eclipse Restart",
    		        "Load Eclipse Restart Files."	,
    		        "Must have created a project prior to loading",
    		        "EclipseLoad20x20.gif");

        StsWorkflowTreeNode defineData = workflowRoot.addChild("Define / Edit", null);
        	StsWorkflowTreeNode horPick = defineData.addChild("com.Sts.Actions.Wizards.Horpick.StsHorpickWizard", "Pick Surfaces",
        		"Auto-track an D event through a seismic volume using the S2S 3D spiral tracker.",
        		"A minimum of one seismic volume is required to pick a surface. Run the Load->3D Seismic Volume first.",
        		"horizonPicker20x20.gif");
        	StsWorkflowTreeNode defineHorizons = defineData.addChild("com.Sts.Actions.Wizards.Horizons.StsHorizonsWizard", "Horizons",
            	"Horizon definition is used to ensure that all surfaces are translated to the same grid and optionally related to markers. It is required for fault definition and reservoir model construction.",
            	"Must have loaded or auto-picked at least one surface. Either run the Load->Surfaces or Define->Pick Surfaces workflow step(s) prior to defining horizons.",
        		"defineHorizons20x20.gif");
//			StsWorkflowTreeNode workflowInitialSeisVel = defineData.addChild("com.Sts.Actions.Wizards.Velocity.StsSeisVelWizard", "Initial Seismic Velocities",
//				"Construct a post-stack velocity model using seismic and well velocity data.",
//				"The minimum reqiurement to construct an initial seismic velocity is either well Time/Depth curve or imported velocity functions.",
//				"loadSeismic20x20.gif");
			StsWorkflowTreeNode workflowVelocityModel = defineData.addChild("com.Sts.Actions.Wizards.Velocity.StsVelocityWizard", "Velocity Model",
        		"Construct a post-stack velocity model using any combination of user input, and imported well markers, horizons and velocity data.",
        		"The minimum reqiurement to construct a velocity model is horizon data or an imported velocity cube. Run the Define->Horizons and/or Load->3D PostStack Volume first. If well ties are desired, import well markers with wells via the Load->Wells & Logs workflow step.",
        		"loadSeismic20x20.gif");
            StsWorkflowTreeNode surfacesFromMarkers = defineData.addChild("com.Sts.Actions.Wizards.SurfacesFromMarkers.StsSurfacesFromMarkersWizard", "Surfaces From Markers",
                 "Construct a surface from a set of markers in existing wells.",
                 "Must have one or more well marker sets in order to construct.",
                 "importSurfaces20x20.gif");
        StsWorkflowTreeNode buildModel = workflowRoot.addChild("Build Model", null);
            StsWorkflowTreeNode defineZones = buildModel.addChild("com.Sts.Actions.Wizards.Zones.StsZonesWizard", "Define Zones",
            	"Select horizon pairs which define zones across the reservoir.",
                "Run the Define->Horizon workflow step first.",
        		"defineHorizons20x20.gif");
            StsWorkflowTreeNode buildBoundary = buildModel.addChild("com.Sts.Actions.Boundary.StsBuildBoundary", "Define Boundary",
            	"Define a boundary within the project to limit the fault and reservoir model construction. Boundaries can be removed and re-defined by running this workflow step.",
            	"Must convert surfaces to horizons prior to defining a boundary. Run the Define->Horizon workflow step first.",
        		"boundary20x20.gif");
        	StsWorkflowTreeNode buildFrame = buildModel.addChild("com.Sts.Actions.Build.StsBuildFrame", "Fault Framework",
            	"Define a fault framework within a pre-defined boundary. Faults are digitized on horizons and/or cursor slices. Dying faults are supportted.",
            	"Must have defined a boundary prior to defining a fault framework. Run the Define->Boundary workflow step first.",
        		"buildFrame20x20.gif");
        	StsWorkflowTreeNode completeModel = buildModel.addChild("com.Sts.Actions.Wizards.Model.StsModelWizard", "Construct Model",
        		"Construct the reservoir model from the fault framework and horizon data and initialize the stratigraphic layering.",
        		"Define boundary, horizons and fault framework prior to constructing the reservoir model. Run Define->Boundary, Define->Horizons and Define->Fault Framework steps in this order prior to constructing the model.",
        		"buildModel20x20.gif");
            StsWorkflowTreeNode exportEclipseModel = buildModel.addChild("com.Sts.Actions.Wizards.SimulationFile.StsSimulationFileWizard", "Export Model in Eclipse Format. ",
        		"Export the reservoir model grid and properties in Eclipse format.",
        		"Define boundary, horizons and fault framework prior to constructing the reservoir model. Run Define->Boundary, Define->Horizons and Define->Fault Framework steps in this order prior to constructing the model.",
        		"buildModel20x20.gif");
            StsWorkflowTreeNode translateEclipse = buildModel.addChild("com.Sts.Actions.Wizards.SimulationFile.StsTranslateEclipseWizard", "Translate Eclipse indexes & oordinates. ",
        		"Translate between S2S IJKB and Eclipse IJK indexes and XYZ coordinates.",
        		"Construct the Eclipse Model first.  Required for index/coordinate translation.",
        		"buildModel20x20.gif");
           StsWorkflowTreeNode eclipseLoad = buildModel.addChild("com.Sts.Actions.Wizards.EclipseLoad.StsEclipseLoadWizard", "Load Eclipse restart file. ",
        		"Reads Eclipse output files and translates to S2S block data structure files, one for each property and restart time (if dynamic).",
        		"Construct the Eclipse Model first.",
        		"buildModel20x20.gif");
/*
        StsWorkflowTreeNode outputData = workflowRoot.addChild("Output", null);
          StsWorkflowTreeNode makeMovie = outputData.addChild("com.Sts.Actions.Wizards.MakeMovie.StsMakeMovieWizard", "Make Movie", "loadAncillaryData20x20.gif");
        treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, makeMovie));
*/
        StsNodeBundle nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { getNode(LOAD_DATA, L_POSTSTACK3D), getNode(PROCESS_SEISMIC, P_POSTSTACK3D) }, StsNodeBundle.ONE_REQUIRED);
        treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, horPick));

        treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, defineHorizons));

        nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { defineHorizons }, StsNodeBundle.ONE_REQUIRED);
        treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, defineZones));
        treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, buildBoundary));

		//nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { getNode(LOAD_DATA, L_POSTSTACK3D), getNode(PROCESS_SEISMIC, P_POSTSTACK3D) }, StsNodeBundle.ONE_REQUIRED);
		//treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, workflowVelocityModel));

        nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { defineHorizons, getNode(LOAD_DATA, L_POSTSTACK3D), getNode(LOAD_DATA, L_SURFACE) }, StsNodeBundle.ONE_REQUIRED);
        treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, workflowVelocityModel));

        nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { getNode(LOAD_DATA, L_WELLSANDLOGS)}, StsNodeBundle.ONE_REQUIRED);
        treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, surfacesFromMarkers));

        nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { buildBoundary }, StsNodeBundle.ONE_REQUIRED);
        treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, buildFrame));
        treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, completeModel));

        nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { completeModel }, StsNodeBundle.ONE_REQUIRED);
        treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, exportEclipseModel));

        nodeBundle = new StsNodeBundle(exportEclipseModel);
        treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, translateEclipse));
        treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, eclipseLoad));
        logUsageChange();
    }

    protected boolean runCreateObjectsPanel(StsObjectTree objectTree, StsModel model)
    {
        try
        {
            if(objectTree == null) return false;
            rootNode = objectTree.createRootNode(model.getProject(), "Project");

            // Data Node
            dataNode = rootNode.addStaticNode("Data");
            seismicNode = checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsSeismicVolume"), "3D Volumes", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsSeismicLineSet"), "2D Line Sets", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsPreStackLineSet3d"), "3D Gathers", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsPreStackLineSet2d"), "2D Gathers", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsPreStackMicroseismicSet"), "Microseismic Gathers", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsWell"), "Wells", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsLogCurveType"), "Log Types", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsTimeLogCurveType"), "Time Log Types", false);

			StsVspClass vspClass = (StsVspClass)model.getCreateStsClass("com.Sts.DBTypes.StsVsp");
            StsTreeNode vspNode = checkAddStaticNode(model, dataNode, vspClass, "VSP", false);
            if(vspNode != null)
            {
			    StsSubType[] subTypes = vspClass.getSubTypes();
			    for(int n = 0; n < subTypes.length; n++)
                    checkAddDynamicNode(model, vspNode, subTypes[n], subTypes[n].getName(), false);
            }
            StsTreeNode virtualVolumeNode = checkAddStaticNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsVirtualVolume"), "Virtual Volumes", false);
            if(virtualVolumeNode != null)
            {
                checkAddDynamicNode(model, virtualVolumeNode, model.getCreateStsClass("com.Sts.DBTypes.StsMathVirtualVolume"), "Math", false);
                checkAddDynamicNode(model, virtualVolumeNode, model.getCreateStsClass("com.Sts.DBTypes.StsBlendedVirtualVolume"), "Blended", false);
                checkAddDynamicNode(model, virtualVolumeNode, model.getCreateStsClass("com.Sts.DBTypes.StsCrossplotVirtualVolume"), "Crossplot", false);
                checkAddDynamicNode(model, virtualVolumeNode, model.getCreateStsClass("com.Sts.DBTypes.StsRGBAVirtualVolume"), "RGBA", false);
                checkAddDynamicNode(model, virtualVolumeNode, model.getCreateStsClass("com.Sts.DBTypes.StsFilterVirtualVolume"), "Filter", false);
                checkAddDynamicNode(model, virtualVolumeNode, model.getCreateStsClass("com.Sts.DBTypes.StsSensorVirtualVolume"), "Sensor", false);
            }
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsCrossplot"), "Crossplot", false);

            subVolumeNode = checkAddStaticNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsSubVolume"), "SubVolumes", false);
            if(subVolumeNode != null)
            {
                checkAddDynamicNode(model, subVolumeNode, model.getCreateStsClass("com.Sts.DBTypes.StsDualSurfaceSubVolume"), "Dual Surface", false);
                checkAddDynamicNode(model, subVolumeNode, model.getCreateStsClass("com.Sts.DBTypes.StsBoxSetSubVolume"), "Box Set", false);
                checkAddDynamicNode(model, subVolumeNode, model.getCreateStsClass("com.Sts.DBTypes.StsWellSubVolume"), "Well Set", false);
            }
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsHorpick"), "Horizon Picks", true);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsSurface"), "Surfaces", true);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsPlatform"), "Drilling Platform", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsWellPlanSet"), "Well Plans", false);

            StsTreeNode sensorNode = checkAddStaticNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsSensor"), "Sensors", false);
            checkAddDynamicNode(model, sensorNode, model.getCreateStsClass("com.Sts.DBTypes.StsStaticSensor"), "Static", false);
            checkAddDynamicNode(model, sensorNode, model.getCreateStsClass("com.Sts.DBTypes.StsDynamicSensor"), "Dynamic", false);

            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsMVFractureSet"), "MV Fractures", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsTriangulatedFracture"), "Fractures", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsMultiAttributeVector"), "Multi-Attribute Vectors", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsFaultStickSet"), "Fault Sticks", true);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsSpectrum"), "Palettes", true);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsCultureObjectSet2D"), "Culture Sets", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsMovie"), "Movies", false);

            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsMonitor"), "Monitors", false);

            StsTreeNode alarmNode = checkAddStaticNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsAlarm"), "Alarm", false);
            checkAddDynamicNode(model, alarmNode, model.getCreateStsClass("com.Sts.DBTypes.StsSurfaceAlarm"), "Surface", false);
            checkAddDynamicNode(model, alarmNode, model.getCreateStsClass("com.Sts.DBTypes.StsWellAlarm"), "Well", false);
            checkAddDynamicNode(model, alarmNode, model.getCreateStsClass("com.Sts.DBTypes.StsValueAlarm"), "Value", false);

            StsAncillaryDataClass adClass = (StsAncillaryDataClass)model.getCreateStsClass("com.Sts.DBTypes.StsAncillaryData");
            StsTreeNode adNode = checkAddStaticNode(model, dataNode, adClass, "Ancillary Data", false);
            if(adNode != null)
            {
			    StsSubType[] subTypes = adClass.getSubTypes();
			    for(int n = 0; n < subTypes.length; n++)
                    checkAddDynamicNode(model, adNode, subTypes[n], subTypes[n].getName(), false);
            }
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsWeighPoint"), "WayPoints", true);


            // Model Node
            modelNode = rootNode.addStaticNode("Model");
            checkAddDynamicNode(model, modelNode, model.getCreateStsClass("com.Sts.DBTypes.StsModelSurface"), "Horizons", true);
            checkAddDynamicNode(model, modelNode, model.getCreateStsClass("com.Sts.DBTypes.StsPreStackVelocityModel3d"), "PreStack3d 3D Velocity Model", false);
            checkAddDynamicNode(model, modelNode, model.getCreateStsClass("com.Sts.DBTypes.StsPreStackVelocityModel2d"), "PreStack3d 2D Velocity Model", false);
            checkAddDynamicNode(model, modelNode, model.getCreateStsClass("com.Sts.DBTypes.StsSeismicVelocityModel"), "Velocity Model", false);
            checkAddDynamicNode(model, modelNode, model.getCreateStsClass("com.Sts.DBTypes.StsEditTdSet"), "Well TD Edits", false);
            checkAddDynamicNode(model, modelNode, model.getCreateStsClass("com.Sts.DBTypes.StsLine"), "Boundary Lines", true);
            checkAddDynamicNode(model, modelNode, model.getCreateStsClass("com.Sts.DBTypes.StsFaultLine"), "Fault Lines", true);
            checkAddDynamicNode(model, modelNode, model.getCreateStsClass("com.Sts.DBTypes.StsZone"), "Units", true);
            checkAddDynamicNode(model, modelNode, model.getCreateStsClass("com.Sts.DBTypes.StsSection"), "Sections", true);
            checkAddDynamicNode(model, modelNode, model.getCreateStsClass("com.Sts.DBTypes.StsFractureSet"), "Fracture Sets", true);
            checkAddDynamicNode(model, modelNode, model.getCreateStsClass("com.Sts.DBTypes.StsBlock"), "Blocks", true);
            checkAddDynamicNode(model, modelNode, model.getCreateStsClass("com.Sts.DBTypes.StsBuiltModel"), "Built Models", true);
            checkAddDynamicNode(model, modelNode, model.getCreateStsClass("com.Sts.DBTypes.StsEclipseModel"), "Eclipse Models", true);
            objectTree.finalizeTreeModel();
            return true;

        }
        catch(Exception e)
        {
            StsException.outputException("StsFractureAnalysisWorkflow.createObjectsPanel() failed.",
                e, StsException.WARNING);
            return false;
        }
}
}
