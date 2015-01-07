package com.Sts.WorkflowPlugIn.PlugIns;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.Types.StsSubType;
import com.Sts.Workflow.*;
import com.Sts.WorkflowPlugIn.*;
import com.Sts.UI.ObjectPanel.StsObjectTree;
import com.Sts.UI.ObjectPanel.StsTreeNode;
import com.Sts.MVC.StsModel;
import com.Sts.DBTypes.StsVspClass;
import com.Sts.DBTypes.StsAncillaryDataClass;
import com.Sts.Utilities.StsException;

public class StsVelocityModelWorkflow extends StsWorkflowPlugIn
{
    public StsVelocityModelWorkflow()
    {
        name = "StsVelocityModelWorkflow";
        workflowName = "Velocity Modeling";
        checkName();
        description = new String("The Velocity Modeling Workflow contains all the steps required" +
                                 " to tie wells and well log information with seismic derived" +
                                 " events such as surfaces. It creates a velocity model" +
                                 " which allows the interactive conversion between time" +
                                 " and depth and is also used for realtime target refinement.");
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
        StsNodeBundle nodeBundle;

        treeModel.addMenuNode(newProject);
        treeModel.addMenuNode(openProject);

        createGroupAndNodes(workflowRoot, treeModel, PROCESS_SEISMIC, new byte[] {P_POSTSTACK3D, P_POSTSTACK2D, P_VSP});
        createGroupAndNodes(workflowRoot, treeModel, LOAD_DATA, new byte[] {L_POSTSTACK3D, L_POSTSTACK2D, L_WELLSANDLOGS, L_VSP, L_SURFACE, L_COLORPALETTE, L_CULTURE, L_ANCILLARY});

        treeModel.addNodeConnection(new StsNodeConnection(getNode(LOAD_DATA, L_WELLSANDLOGS), getNode(LOAD_DATA, L_VSP)));

        StsWorkflowTreeNode defineData = workflowRoot.addChild("Define / Edit", null);
    		StsWorkflowTreeNode workflowEditWellMarkers = defineData.addChild("com.Sts.Actions.Wizards.EditWellMarkers.StsEditWellMarkersWizard", "Well Markers",
              	"Edit previously loaded well markers or interpret new markers on previously loaded wells and logs. Markers are automatically loaded with associated wells when the wells are loaded.",
              	"Wells, logs and possibly existing markers must have already been loaded prior to editing well markers. Run the Load->Wells & Logs workflow step first.",
    			"markersOnSurface20x20.gif");
       	    StsWorkflowTreeNode workflowEditLogCurves = defineData.addChild("com.Sts.Actions.Wizards.LogEdit.StsEditLogCurveWizard", "Log Curves",
        		"Create and/or edit log curves associated with each well.",
        		" Run Load->Wells & Logs first.",
        		"tdEdit20x20.gif");
            StsWorkflowTreeNode workflowEditTdCurves = defineData.addChild("com.Sts.Actions.Wizards.EditTd.StsEditTdCurveWizard", "Time-Depth Curves",
        		"Create and/or edit time-depth curves associated with each well. The well is displayed alongside any seismic attribute data and markers can be re-positioned relative to the attribute data, effectively adjusting the time-depth function.",
        		"Seismic attribute data and all well, marker and existing time-depth functions must be loaded prior to editting time-depth curves. Run Load->Wells & Logs and Load->3D PostStack Seismic first.",
        		"tdEdit20x20.gif");
        	StsWorkflowTreeNode workflowDefineHorizons = defineData.addChild("com.Sts.Actions.Wizards.Horizons.StsHorizonsWizard", "Horizons",
        		"Horizon definition is used to ensure that all surfaces are translated to the same grid and optionally related to markers. It is required for fault definition and reservoir model construction.",
        		"Must have loaded or auto-picked at least one surface. Either run the Load->Surfaces or Define->Pick Surfaces workflow step(s) prior to defining horizons.",
        		"defineHorizons20x20.gif");
			StsWorkflowTreeNode workflowInitialSeisVel = defineData.addChild("com.Sts.Actions.Wizards.Velocity.StsSeisVelWizard", "Initial Seismic Velocities",
				"Construct a post-stack velocity model using seismic and well velocity data.",
				"The minimum reqiurement to construct an initial seismic velocity is either well Time/Depth curve or imported velocity functions.",
				"loadSeismic20x20.gif");
			StsWorkflowTreeNode workflowVelocityModel = defineData.addChild("com.Sts.Actions.Wizards.Velocity.StsVelocityWizard", "Velocity Model",
            	"Construct a post-stack velocity model using any combination of user input, and imported well markers, horizons and velocity data.",
            	"The minimum reqiurement to construct a velocity model is horizon data or an imported velocity cube. Run the Define->Horizons and/or Load->3D PostStack Volume first. If well ties are desired, import well markers with wells via the Load->Wells & Logs workflow step.",
        		"loadSeismic20x20.gif");
        	StsWorkflowTreeNode filterVolume = defineData.addChild("com.Sts.Actions.Wizards.VirtualVolume.StsFilterVolumeWizard", "Filter Volume",
        		"Filter an existing post-stack volume. Filters primarily consist of smoothing operators with user defined kernal types and sizes. Filter volumes are a type of virtual volume.",
            	"Must have a post-stack seismic volume in order to run the filter volume workflow step. Poststack volume can be directly loaded or generated by building a velocity modeling.",
        		"defineVirtual20x20.gif");
            StsWorkflowTreeNode surfacesFromMarkers = defineData.addChild("com.Sts.Actions.Wizards.SurfacesFromMarkers.StsSurfacesFromMarkersWizard", "Surfaces From Markers",
                 "Construct a surface from a set of markers in existing wells.",
                 "Must have one or more well marker sets in order to construct.",
                 "importSurfaces20x20.gif");

            nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { getNode(LOAD_DATA, L_WELLSANDLOGS), surfacesFromMarkers }, StsNodeBundle.ONE_REQUIRED);
            treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, workflowEditWellMarkers));

        	nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { getNode(LOAD_DATA, L_WELLSANDLOGS), surfacesFromMarkers }, StsNodeBundle.ONE_REQUIRED);
        	treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, workflowEditTdCurves));

            treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, workflowDefineHorizons));

	        nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { getNode(LOAD_DATA, L_POSTSTACK3D), getNode(PROCESS_SEISMIC, P_POSTSTACK3D)}, StsNodeBundle.ONE_REQUIRED);
	        treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, workflowInitialSeisVel));

            nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { workflowInitialSeisVel, getNode(LOAD_DATA, L_POSTSTACK3D), workflowDefineHorizons}, StsNodeBundle.ONE_REQUIRED);
            treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, workflowVelocityModel));

			nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { workflowVelocityModel, getNode(LOAD_DATA, L_POSTSTACK3D), getNode(LOAD_DATA, L_SURFACE)}, StsNodeBundle.ONE_REQUIRED);
			treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, filterVolume));

			nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { getNode(LOAD_DATA, L_WELLSANDLOGS)}, StsNodeBundle.ONE_REQUIRED);
            treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, surfacesFromMarkers));
/*
        StsWorkflowTreeNode outputData = workflowRoot.addChild("Output", null);
          	StsWorkflowTreeNode makeMovie = outputData.addChild("com.Sts.Actions.Wizards.MakeMovie.StsMakeMovieWizard", "QT Movie", "loadAncillaryData20x20.gif");
          	treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, makeMovie));
*/
        logUsageChange();
    }

    public void addOptionalNodes(StsTreeModel treeModel, String[] options)
    {

    }

    protected boolean runCreateObjectsPanel(StsObjectTree objectTree, StsModel model)
    {
        try
        {
            if(objectTree == null) return false;
            rootNode = objectTree.createRootNode(model.getProject(), "Project");

            // Data Node
            dataNode = rootNode.addStaticNode("Data");
            seismicNode = checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsSeismicVolume"), "3D Volumes", true);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsSeismicLineSet"), "2D Line Sets", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsPreStackLineSet3d"), "3D Gathers", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsPreStackLineSet2d"), "2D Gathers", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsPreStackMicroseismicSet"), "Microseismic Gathers", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsWell"), "Wells", true);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsLiveWell"), "Live Wells", false);            
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsLogCurveType"), "Log Types", true);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsTimeLogCurveType"), "Time Log Types", false);

			StsVspClass vspClass = (StsVspClass)model.getCreateStsClass("com.Sts.DBTypes.StsVsp");
            StsTreeNode vspNode = checkAddStaticNode(model, dataNode, vspClass, "VSP", true);
            if(vspNode != null)
            {
			    StsSubType[] subTypes = vspClass.getSubTypes();
			    for(int n = 0; n < subTypes.length; n++)
                    checkAddDynamicNode(model, vspNode, subTypes[n], subTypes[n].getName(), true);
            }
            StsTreeNode virtualVolumeNode = checkAddStaticNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsVirtualVolume"), "Virtual Volumes", true);
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

            subVolumeNode = checkAddStaticNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsSubVolume"), "SubVolumes", true);
            if(subVolumeNode != null)
            {
                checkAddDynamicNode(model, subVolumeNode, model.getCreateStsClass("com.Sts.DBTypes.StsDualSurfaceSubVolume"), "Dual Surface", false);
                checkAddDynamicNode(model, subVolumeNode, model.getCreateStsClass("com.Sts.DBTypes.StsBoxSetSubVolume"), "Box Set", false);
                checkAddDynamicNode(model, subVolumeNode, model.getCreateStsClass("com.Sts.DBTypes.StsWellSubVolume"), "Well Set", false);
            }
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsHorpick"), "Horizon Picks", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsSurface"), "Surfaces", true);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsPlatform"), "Drilling Platform", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsWellPlanSet"), "Well Plans", false);

            StsTreeNode sensorNode = checkAddStaticNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsSensor"), "Sensors", true);
            checkAddDynamicNode(model, sensorNode, model.getCreateStsClass("com.Sts.DBTypes.StsStaticSensor"), "Static", false);
            checkAddDynamicNode(model, sensorNode, model.getCreateStsClass("com.Sts.DBTypes.StsDynamicSensor"), "Dynamic", false);

            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsMVFractureSet"), "MV Fractures", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsTriangulatedFracture"), "Fractures", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsMultiAttributeVector"), "Multi-Attribute Vectors", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsFaultStickSet"), "Fault Sticks", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsSpectrum"), "Palettes", true);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsCultureObjectSet2D"), "Culture Sets", true);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsMovie"), "Movies", false);

            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsMonitor"), "Monitors", false);

            StsTreeNode alarmNode = checkAddStaticNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsAlarm"), "Alarm", true);
            checkAddDynamicNode(model, alarmNode, model.getCreateStsClass("com.Sts.DBTypes.StsSurfaceAlarm"), "Surface", false);
            checkAddDynamicNode(model, alarmNode, model.getCreateStsClass("com.Sts.DBTypes.StsWellAlarm"), "Well", false);
            checkAddDynamicNode(model, alarmNode, model.getCreateStsClass("com.Sts.DBTypes.StsValueAlarm"), "Value", false);

            StsAncillaryDataClass adClass = (StsAncillaryDataClass)model.getCreateStsClass("com.Sts.DBTypes.StsAncillaryData");
            StsTreeNode adNode = checkAddStaticNode(model, dataNode, adClass, "Ancillary Data", true);
            if(adNode != null)
            {
			    StsSubType[] subTypes = adClass.getSubTypes();
			    for(int n = 0; n < subTypes.length; n++)
                    checkAddDynamicNode(model, adNode, subTypes[n], subTypes[n].getName(), true);
            }
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsWeighPoint"), "WayPoints", true);


            // Model Node
            modelNode = rootNode.addStaticNode("Model");
            checkAddDynamicNode(model, modelNode, model.getCreateStsClass("com.Sts.DBTypes.StsModelSurface"), "Horizons", true);
            checkAddDynamicNode(model, modelNode, model.getCreateStsClass("com.Sts.DBTypes.StsPreStackVelocityModel3d"), "PreStack3d 3D Velocity Model", true);
            checkAddDynamicNode(model, modelNode, model.getCreateStsClass("com.Sts.DBTypes.StsPreStackVelocityModel2d"), "PreStack3d 2D Velocity Model", true);
            checkAddDynamicNode(model, modelNode, model.getCreateStsClass("com.Sts.DBTypes.StsSeismicVelocityModel"), "Velocity Model", true);
            checkAddDynamicNode(model, modelNode, model.getCreateStsClass("com.Sts.DBTypes.StsEditTdSet"), "Well TD Edits", true);
            checkAddDynamicNode(model, modelNode, model.getCreateStsClass("com.Sts.DBTypes.StsLine"), "Boundary Lines", true);
            checkAddDynamicNode(model, modelNode, model.getCreateStsClass("com.Sts.DBTypes.StsFaultLine"), "Fault Lines", true);
            checkAddDynamicNode(model, modelNode, model.getCreateStsClass("com.Sts.DBTypes.StsZone"), "Units", true);
            checkAddDynamicNode(model, modelNode, model.getCreateStsClass("com.Sts.DBTypes.StsSection"), "Sections", true);
            checkAddDynamicNode(model, modelNode, model.getCreateStsClass("com.Sts.DBTypes.StsFractureSet"), "Fracture Sets", true);
            checkAddDynamicNode(model, modelNode, model.getCreateStsClass("com.Sts.DBTypes.StsBlock"), "Blocks", true);

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
