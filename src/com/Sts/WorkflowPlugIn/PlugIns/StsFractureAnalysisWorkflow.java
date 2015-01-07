package com.Sts.WorkflowPlugIn.PlugIns;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.MVC.View3d.*;
import com.Sts.MVC.StsModel;
import com.Sts.Workflow.*;
import com.Sts.WorkflowPlugIn.*;
import com.Sts.UI.ObjectPanel.StsObjectTree;
import com.Sts.UI.ObjectPanel.StsTreeNode;
import com.Sts.DBTypes.StsVspClass;
import com.Sts.DBTypes.StsAncillaryDataClass;
import com.Sts.Types.StsSubType;
import com.Sts.Utilities.StsException;

public class StsFractureAnalysisWorkflow extends StsWorkflowPlugIn
{
    public StsFractureAnalysisWorkflow()
    {
        name = "StsFractureAnalysisWorkflow";
        workflowName = "Fracture Analysis";
        checkName();
        description = new String("The Fracture Analysis Workflow contains all the steps required to load wells," +
                                 " logs, seismic attributes, and sensor data to analyze natural and hydraulic" +
                                 " fracture systems. Also, included is the ability to analyze hydrualic fracture results" +
                                 " in real-time with full control to set current model time, time ranges and playback" +
                                 " speeds, resulting in the most complete understanding of fracture dynamics." );
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

        createGroupAndNodes(workflowRoot, treeModel, PROCESS_SEISMIC, new byte[] {P_POSTSTACK2D, P_POSTSTACK3D, P_VSP});
        createGroupAndNodes(workflowRoot, treeModel, LOAD_DATA, new byte[] {L_POSTSTACK2D, L_POSTSTACK3D, L_WELLSANDLOGS, L_VSP, L_SENSOR, L_SURFACE, L_COLORPALETTE, L_CULTURE, L_ANCILLARY});
		StsWorkflowTreeNode loadMVFractures = getGroupNode(LOAD_DATA).addChild("com.Sts.Actions.Wizards.FracSetMVLoad.StsMVFracSetLoadWizard", "MV Fractures",
                    "Load Fracture Sets created in Midland Valley's 4DFrac",
                    "Must have created an S2S Project.",
                    "MVFractures20x20.gif");
		StsWorkflowTreeNode loadGolderFractures = getGroupNode(LOAD_DATA).addChild("com.Sts.Actions.Wizards.FracSetGolderLoad.StsGolderFracSetLoadWizard", "Golder Fractures",
                    "Load Fracture Sets created in Golder's Fracman",
                    "Must have created an S2S Project.",
                    "MVFractures20x20.gif");
        StsWorkflowTreeNode faultSticks = getGroupNode(LOAD_DATA).addChild("com.Sts.Actions.Wizards.FaultSticks.StsFaultSticksWizard", "Fault Sticks",
	    		"Load ASCII file containing fault sticks"	,
	    		"Must have created an S2S Project.",
	    		"faultSticks20x20.gif");
    	StsWorkflowTreeNode loadStsObjects = getGroupNode(LOAD_DATA).addChild("com.Sts.Actions.Wizards.StsObjectLoad.StsObjectLoadWizard", "S2S Objects",
            	"Load S2S Objects that were exported from other S2S Projects",
            	"Must have created an S2S Project.",
            	"Proximity20x20.gif");

        StsWorkflowTreeNode defineGroup = workflowRoot.addChild("Define / Edit", null);
    	    StsWorkflowTreeNode subVolume = defineGroup.addChild("com.Sts.Actions.Wizards.SubVolume.StsSubVolumeWizard", "SubVolumes",
        		"Define a sub-volume of the current project. Sub-volumes can be defined in relation to surfaces or as an arbitrary set of polyhedron. Sub-volumes can be applied to any seismic volume or crossplot to limit the amount of data displayed on the cursor planes or crossplotted."		,
        		"A minimum of one seismic volume must exist prior to defining a subvolume. Run the Load->3D Poststack Seismic before defining a subvolume.",
        		"defineCropVolume20x20.gif");
            StsWorkflowTreeNode workflowDefineHorizons = defineGroup.addChild("com.Sts.Actions.Wizards.Horizons.StsHorizonsWizard", "Horizons",
        		"Horizon definition is used to ensure that all surfaces are translated to the same grid and optionally related to markers. It is required for fault definition and reservoir model construction.",
        		"Must have loaded or auto-picked at least one surface. Either run the Load->Surfaces or Define->Pick Surfaces workflow step(s) prior to defining horizons.",
        		"defineHorizons20x20.gif");
           StsWorkflowTreeNode griddedSensorVolume = defineGroup.addChild("com.Sts.Actions.Wizards.GriddedSensorAtts.StsGriddedSensorAttsWizard", "Gridded Sensor Attributes",
        		"Create a gridded sensor attribute volume from sensor event data. These volumes are computed on the fly from sensor event data. Virtual volumes are treated as real volumes everywhere else in S2S."	,
        		"Must have at least one sensor event set loaded to define a gridded sensor volume. Run the Load->Sensor workflow step first.",
        		"defineVirtual20x20.gif");
           StsWorkflowTreeNode sensorPartitions = defineGroup.addChild("com.Sts.Actions.Wizards.SensorPartition.StsSensorPartitionWizard", "Partition a Sensor",
        		"Break a single sensor into multiple sensors based on time increments."	,
        		"Must have at least one sensor event set loaded. Run the Load->Sensor workflow step first.",
        		"sensorPartition.gif");
			StsWorkflowTreeNode workflowEditWellMarkers = defineGroup.addChild("com.Sts.Actions.Wizards.EditWellMarkers.StsEditWellMarkersWizard", "Well Markers",
              	"Edit previously loaded well markers or interpret new markers on previously loaded wells and logs. Markers are automatically loaded with associated wells when the wells are loaded.",
              	"Wells, logs and possibly existing markers must have already been loaded prior to editing well markers. Run the Load->Wells & Logs workflow step first.",
    			"markersOnSurface20x20.gif");
            StsWorkflowTreeNode fractures = defineGroup.addChild("com.Sts.Actions.Wizards.Fracture.StsFractureWizard", "Fractures",
        		"Interpret fractures by digitizing them in 3D on the Z slice. Guide the picking with seismic, sensors and surfaces.",
            	"Must have a 3D seismic volume prior to interpreting a fracture set. Run the Load->3D PostStack Seismic workflow step first.",
          		"edgeOnCursor.gif");
        	nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { getNode(LOAD_DATA, L_POSTSTACK3D), getNode(PROCESS_SEISMIC, P_POSTSTACK3D) }, StsNodeBundle.ONE_REQUIRED);
            treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, fractures));
        	nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { getNode(LOAD_DATA, L_WELLSANDLOGS)}, StsNodeBundle.ONE_REQUIRED);
            treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, workflowEditWellMarkers));
        	treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, griddedSensorVolume));
            treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, sensorPartitions));
        	treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, faultSticks));
            treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, loadMVFractures));
			treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, loadGolderFractures));
        	treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, loadStsObjects));
            treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, subVolume));
            treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, workflowDefineHorizons));

        StsWorkflowTreeNode analysisGroup = workflowRoot.addChild("Analyze", null);
        	StsWorkflowTreeNode multiAttributeAnalysis = analysisGroup.addChild("com.Sts.Actions.Wizards.MultiAttributeAnalysis.StsMultiAttrAnalysisWizard", "3D Attribute Vector",
        		"Create a dynamic map of vectors on the Z slice or on a surface, where the azimuth, length and color are controlled by three seismic volumes. This capability is often used to analyze natural fracture patterns.",
        		"Must have loaded or created at least two seismic volumes to run the 3D Vector Workflow Step. Run Load->3D PostStack Seismic or Define->Seismic Attributes (Seismic Attribute Workflow) prior to running 3D Attribute Vector",
        		"3DVectors20x20.gif");
        	nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { getNode(LOAD_DATA, L_POSTSTACK3D) }, StsNodeBundle.ONE_REQUIRED);
        	treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, multiAttributeAnalysis));
        	StsWorkflowTreeNode clusterAnalysis = analysisGroup.addChild("com.Sts.Actions.Wizards.ClusterAnalysis.StsClusterAnalysisWizard", "Cluster Analysis",
            	"Analyze sensor events relative to each other in time and space.",
            	"Must have loaded at least one sensor file",
            	"Cluster20x20.gif");
        	StsWorkflowTreeNode proximityAnalysis = analysisGroup.addChild("com.Sts.Actions.Wizards.ProximityAnalysis.StsProximityAnalysisWizard", "Proximity Analysis",
                	"Analyze sensor event distances relative to a well.",
                	"Must have loaded at least one sensor file and one well",
                	"Proximity20x20.gif");
        	StsWorkflowTreeNode sensorCompare = analysisGroup.addChild("com.Sts.Actions.Wizards.SensorCompare.StsSensorCompareWizard", "Sensor Comparison",
                	"Compare events collected from multiple sensor arrays and configurations.",
                	"Must have loaded at least two sensor datasets. This step compares the events between the sensor datasets.",
                	"SensorCompare20x20.gif");
        	StsWorkflowTreeNode sensorXplot = analysisGroup.addChild("com.Sts.Actions.Wizards.SensorXplot.StsSensorXplotWizard", "Sensor Crossplot",
                	"Conduct crossplot analysis of sensor attributes.",
                	"Must have loaded at least one sensor dataset. This step allows the comparison of sensor attributes by plotting those attributes against one another.",
                	"Xplot20x20.gif");
        	StsWorkflowTreeNode sensorQualify = analysisGroup.addChild("com.Sts.Actions.Wizards.SensorQualify.StsSensorQualifyWizard", "Sensor Qualification",
                	"Limit sensor events to specific attribute ranges.",
                	"Must have loaded at least one sensor dataset. This step allows the qualification of sensor events for further analysis.",
                	"sensorQualify.gif");
        	StsWorkflowTreeNode volumeStimulated = analysisGroup.addChild("com.Sts.Actions.Wizards.VolumeStimulated.StsVolumeStimulatedWizard", "Volume Stimulated",
                	"Determine the volume of reservoir stimulated.",
                	"Must have loaded at least one sensor dataset.",
                	"Volumetrics20x20.gif");
         	StsWorkflowTreeNode surfaceCurvature = analysisGroup.addChild("com.Sts.Actions.Wizards.SurfaceCurvature.StsSurfaceCurvatureWizard", "Surface Curvature",
                	"Compute Curvature Attributes on a Surface",
                	"Must have loaded at least one surface.",
                	"importSurfaces20x20.gif");
            StsWorkflowTreeNode perforationAssignment = analysisGroup.addChild("com.Sts.Actions.Wizards.PerforationAttributes.StsPerforationAttributesWizard", "Perforation Assignment",
                	"Assign Perforations to Sensors and Compute Attributes",
                	"One Well and One Sensor",
                	"sensorQualify.gif");
        	StsWorkflowTreeNode sensorCorrelation = analysisGroup.addChild("com.Sts.Actions.Wizards.SensorCorrelation.StsSensorCorrelationWizard", "Sensor Correlation",
                	"Correlate sensor with other sensors",
                	"Must have loaded at least two sensor datasets. This step allows the correlation of sensor events with other related sensor data.",
                	"sensorCorrelation20x20.gif");

            nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { getNode(LOAD_DATA, L_SENSOR) }, StsNodeBundle.ONE_REQUIRED);
            treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, sensorQualify));
        	treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, sensorCompare));
            treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, sensorXplot));
        	treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, clusterAnalysis));
        	treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, volumeStimulated));
        	treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, proximityAnalysis));
            treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, perforationAssignment));
            treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, sensorCorrelation));
            treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, surfaceCurvature));

        StsWorkflowTreeNode outputData = workflowRoot.addChild("Output", null);
            StsWorkflowTreeNode makeMovie = outputData.addChild("com.Sts.Actions.Wizards.MakeMovie.StsMakeMovieWizard", "QT Movie", "loadAncillaryData20x20.gif");
            treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, makeMovie));

        logUsageChange();
    }

    public void addAdditionalToolbars(StsWin3dBase win3d)
    {
        //if(win3d instanceof StsWin3d)
        //    win3d.checkAddTimeActionToolbar();
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
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsSeismicLineSet"), "2D Line Sets", true);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsPreStackLineSet3d"), "3D Gathers", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsPreStackLineSet2d"), "2D Gathers", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsPreStackMicroseismicSet"), "Microseismic Gathers", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsLiveWell"), "Live Wells", false);                        
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsWell"), "Wells", true);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsLogCurveType"), "Log Types", true);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsTimeLogCurveType"), "Time Log Types", false);

			StsVspClass vspClass = (StsVspClass)model.getCreateStsClass("com.Sts.DBTypes.StsVsp");
            StsTreeNode vspNode = checkAddStaticNode(model, dataNode, vspClass, "VSP", false);
            if(vspNode != null)
            {
			    StsSubType[] subTypes = vspClass.getSubTypes();
			    for(int n = 0; n < subTypes.length; n++)
                    checkAddDynamicNode(model, vspNode, subTypes[n], subTypes[n].getName(), false);
            }
            StsTreeNode virtualVolumeNode = checkAddStaticNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsVirtualVolume"), "Virtual Volumes", true);
            if(virtualVolumeNode != null)
            {
                checkAddDynamicNode(model, virtualVolumeNode, model.getCreateStsClass("com.Sts.DBTypes.StsMathVirtualVolume"), "Math", false);
                checkAddDynamicNode(model, virtualVolumeNode, model.getCreateStsClass("com.Sts.DBTypes.StsBlendedVirtualVolume"), "Blended", false);
                checkAddDynamicNode(model, virtualVolumeNode, model.getCreateStsClass("com.Sts.DBTypes.StsCrossplotVirtualVolume"), "Crossplot", false);
                checkAddDynamicNode(model, virtualVolumeNode, model.getCreateStsClass("com.Sts.DBTypes.StsRGBAVirtualVolume"), "RGBA", false);
                checkAddDynamicNode(model, virtualVolumeNode, model.getCreateStsClass("com.Sts.DBTypes.StsFilterVirtualVolume"), "Filter", false);
                checkAddDynamicNode(model, virtualVolumeNode, model.getCreateStsClass("com.Sts.DBTypes.StsSensorVirtualVolume"), "Sensor", true);
            }
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsCrossplot"), "Crossplot", false);

            subVolumeNode = checkAddStaticNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsSubVolume"), "SubVolumes", true);
            if(subVolumeNode != null)
            {
                checkAddDynamicNode(model, subVolumeNode, model.getCreateStsClass("com.Sts.DBTypes.StsDualSurfaceSubVolume"), "Dual Surface", true);
                checkAddDynamicNode(model, subVolumeNode, model.getCreateStsClass("com.Sts.DBTypes.StsBoxSetSubVolume"), "Box Set", false);
                checkAddDynamicNode(model, subVolumeNode, model.getCreateStsClass("com.Sts.DBTypes.StsWellSubVolume"), "Well Set", false);
            }
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsHorpick"), "Horizon Picks", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsSurface"), "Surfaces", true);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsPlatform"), "Drilling Platform", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsWellPlanSet"), "Well Plans", false);

            StsTreeNode sensorFilterNode = checkAddStaticNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsSensorFilter"), "Sensor Filters", true);            
                checkAddDynamicNode(model, sensorFilterNode, model.getCreateStsClass("com.Sts.DBTypes.StsSensorQualifyFilter"), "Qualify", true);
                checkAddDynamicNode(model, sensorFilterNode, model.getCreateStsClass("com.Sts.DBTypes.StsSensorProximityFilter"), "Proximity", true);
                checkAddDynamicNode(model, sensorFilterNode, model.getCreateStsClass("com.Sts.DBTypes.StsSensorPartitionFilter"), "Partition", true);

            StsTreeNode sensorNode = checkAddStaticNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsSensor"), "Sensors", true);            
                checkAddDynamicNode(model, sensorNode, model.getCreateStsClass("com.Sts.DBTypes.StsStaticSensor"), "Static", true);
                checkAddDynamicNode(model, sensorNode, model.getCreateStsClass("com.Sts.DBTypes.StsDynamicSensor"), "Dynamic", true);

            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsMVFractureSet"), "MV Fractures", true);
			checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsGolderFractureSet"), "Golder Fractures", true);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsTriangulatedFracture"), "Fractures", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsMultiAttributeVector"), "Multi-Attribute Vectors", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsFaultStickSet"), "Fault Sticks", true);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsSpectrum"), "Palettes", true);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsCultureObjectSet2D"), "Culture Sets", true);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsMovie"), "Movies", true);

            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsMonitor"), "Monitors", false);

            StsTreeNode alarmNode = checkAddStaticNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsAlarm"), "Alarm", false);
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
