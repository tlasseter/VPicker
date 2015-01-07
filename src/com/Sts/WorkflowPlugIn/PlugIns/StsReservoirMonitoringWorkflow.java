package com.Sts.WorkflowPlugIn.PlugIns;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.Utilities.*;
import com.Sts.Workflow.*;
import com.Sts.WorkflowPlugIn.*;

public class StsReservoirMonitoringWorkflow extends StsWorkflowPlugIn
{
    public StsReservoirMonitoringWorkflow()
    {
        name = "StsReservoirMonitoringWorkflow";
        workflowName = "Reservoir Monitoring";
        checkName();
        description = new String("The Reservoir Monitoring Workflow contains all the steps required to load wells," +
                                 " logs, seismic attributes, reservoir models, reservoir properties and simulation" +
                                 " results into the system and correlate real-time data with them. Also, included is" +
                                 " all the velocity modeling capabilities required to allow the depth conversion" +
                                 " of any time or seismic depth data." );
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

        createGroupAndNodes(workflowRoot, treeModel, LOAD_DATA, new byte[] {L_WELLSANDLOGS, L_SENSOR, L_SURFACE });


		    StsWorkflowTreeNode loadDtsRates = getGroupNode(LOAD_DATA).addChild("com.Sts.Actions.Wizards.DTS.BhDtsLoadWizard", "Temperature",
    		        "Load Distributed Temperature Sensor Files."	,
    		        "Must have loaded wells to associate with the Temperature Sensor data",
    		        "TempSensors20x20.gif");

        StsWorkflowTreeNode defineGroup = workflowRoot.addChild("Define / Edit", null);
           StsWorkflowTreeNode workflowMonitor = defineGroup.addChild("com.Sts.Actions.Wizards.Monitor.StsMonitorWizard", "Monitors",
        		"Define a data monitor for dynamic data collection. Data monitors are objects that watch directories, files and connections for new data and automatically load them into the database.",
        	   	"Need to open or create a project prior to defining a data monitor.",
          		"monitor20x20.gif");

           StsWorkflowTreeNode alarms = defineGroup.addChild("com.Sts.Actions.Wizards.Alarms.StsAlarmsWizard", "Alarms",
        		"Define alarms for realtime data loaders.",
        	   	"Need to have loaded sensor data prior to defining alarms.",
          		"alarms20x20.gif");

           StsWorkflowTreeNode movie = defineGroup.addChild("com.Sts.Actions.Wizards.Movie.StsMovieWizard", "Animations",
          	    "Animate the cursor planes and user views. control the speed and interval a slice is moved through the volume while rotating the users azimuth and elevation relative to the data.",
          	    "No pre-requisites to create and run an animation. May not be very interesting without data, but it will work.",
          	  	"createMovie20x20.gif");

           StsWorkflowTreeNode surfacesFromMarkers = defineGroup.addChild("com.Sts.Actions.Wizards.SurfacesFromMarkers.StsSurfacesFromMarkersWizard", "Surfaces From Markers",
                 "Construct a surface from a set of markers in existing wells.",
                 "Must have one or more well marker sets in order to construct.",
                 "importSurfaces20x20.gif");

           StsWorkflowTreeNode selectActiveWellZone = defineGroup.addChild("com.Sts.Actions.Wizards.ActiveWell.StsActiveWellWizard", "Active Well Zone",
                 "Select the active well zone and highlight in views",
                 "Must have one well",
                 "activeWells20x20.gif");

           treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, workflowMonitor));
           treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, alarms));
           treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, movie));
           nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { getNode(LOAD_DATA, L_WELLSANDLOGS) }, StsNodeBundle.ONE_REQUIRED);
           treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, loadDtsRates));
           nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { getNode(LOAD_DATA, L_WELLSANDLOGS), surfacesFromMarkers }, StsNodeBundle.ONE_REQUIRED);
           treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, surfacesFromMarkers));
           treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, selectActiveWellZone));

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
            seismicNode = checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsSeismicVolume"), "3D Volumes", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsSeismicLineSet"), "2D Line Sets", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsPreStackLineSet3d"), "3D Gathers", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsPreStackLineSet2d"), "2D Gathers", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsPreStackMicroseismicSet"), "Microseismic Gathers", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsWell"), "Wells", true);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsLiveWell"), "Live Wells", true);            
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsLogCurveType"), "Log Types", true);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsTimeLogCurveType"), "Time Log Types", true);
            //checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsTank"), "Storage Tanks", false);
            //checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsPump"), "Pumps", false);

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

            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsMVFractureSet"), "MV Fractures", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsTriangulatedFracture"), "Fractures", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsMultiAttributeVector"), "Multi-Attribute Vectors", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsFaultStickSet"), "Fault Sticks", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsSpectrum"), "Palettes", true);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsCultureObjectSet2D"), "Culture Sets", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsMovie"), "Movies", true);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsMonitor"), "Monitors", true);

            StsTreeNode alarmNode = checkAddStaticNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsAlarm"), "Alarms", true);
            checkAddDynamicNode(model, alarmNode, model.getCreateStsClass("com.Sts.DBTypes.StsSurfaceAlarm"), "Surface", true);
            checkAddDynamicNode(model, alarmNode, model.getCreateStsClass("com.Sts.DBTypes.StsWellAlarm"), "Well", true);
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
            checkAddDynamicNode(model, modelNode, model.getCreateStsClass("com.Sts.DBTypes.StsModelSurface"), "Horizons", false);
            checkAddDynamicNode(model, modelNode, model.getCreateStsClass("com.Sts.DBTypes.StsPreStackVelocityModel3d"), "PreStack3d 3D Velocity Model", false);
            checkAddDynamicNode(model, modelNode, model.getCreateStsClass("com.Sts.DBTypes.StsPreStackVelocityModel2d"), "PreStack3d 2D Velocity Model", false);
            checkAddDynamicNode(model, modelNode, model.getCreateStsClass("com.Sts.DBTypes.StsSeismicVelocityModel"), "Velocity Model", false);
            checkAddDynamicNode(model, modelNode, model.getCreateStsClass("com.Sts.DBTypes.StsEditTdSet"), "Well TD Edits", false);
            checkAddDynamicNode(model, modelNode, model.getCreateStsClass("com.Sts.DBTypes.StsLine"), "Boundary Lines", false);
            checkAddDynamicNode(model, modelNode, model.getCreateStsClass("com.Sts.DBTypes.StsFaultLine"), "Fault Lines", false);
            checkAddDynamicNode(model, modelNode, model.getCreateStsClass("com.Sts.DBTypes.StsZone"), "Units", false);
            checkAddDynamicNode(model, modelNode, model.getCreateStsClass("com.Sts.DBTypes.StsSection"), "Sections", false);
            checkAddDynamicNode(model, modelNode, model.getCreateStsClass("com.Sts.DBTypes.StsFractureSet"), "Fracture Sets", false);
            checkAddDynamicNode(model, modelNode, model.getCreateStsClass("com.Sts.DBTypes.StsBlock"), "Blocks", false);            

            objectTree.finalizeTreeModel();
            return true;

        }
        catch(Exception e)
        {
            StsException.outputException("StsReservoirMonitoringWorkflow.createObjectsPanel() failed.",
                e, StsException.WARNING);
            return false;
        }
    }
}
