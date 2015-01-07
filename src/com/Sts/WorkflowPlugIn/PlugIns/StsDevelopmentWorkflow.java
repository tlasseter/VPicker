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
import com.Sts.Workflow.*;
import com.Sts.WorkflowPlugIn.*;

public class StsDevelopmentWorkflow extends StsWorkflowPlugIn
{
    public StsDevelopmentWorkflow()
    {
        name = "StsDevelopmentWorkflow";
        workflowName = "Development Workflow";
        checkName();
        description = new String("The Development Workflow contains workflow steps that are under development" +
        		"and additionally steps that may be required to execute the steps under development.");
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

        StsWorkflowTreeNode loadGroup = workflowRoot.addChild("Load Data", null);
		    StsWorkflowTreeNode openSpiritWells = loadGroup.addChild("com.Sts.Actions.Wizards.OSWell.StsOSWellWizard", "Wells",
    		        "Load well and well related data from an OpenSpirit datastore."	,
    		        "Must have established an OpenSpirit connection prior to attempting to load data.",
    		        "openSpiritWells20x20.gif");

		    StsWorkflowTreeNode loadMicroseismicPrestack = loadGroup.addChild("com.Sts.Actions.Wizards.MicroseismicPreStack.StsMicroPreStackWizard", "Microseismic Gathers",
	    		    "Load raw microseismic data from SegY files.",
	    		    "Must have created a project. Once loaded must have related events to access gathers.",
	    		    "2DSeismic20x20.gif");

		    StsWorkflowTreeNode loadDtsRates = loadGroup.addChild("com.Sts.Actions.Wizards.DTS.BhDtsLoadWizard", "Temperature",
    		        "Load Distributed Temperature Sensor Files."	,
    		        "Must have loaded wells to associate with the Temperature Sensor data",
    		        "TempSensors20x20.gif");

		    StsWorkflowTreeNode loadPumpRates = loadGroup.addChild("com.Sts.Actions.Wizards.FlowSystem.PumpRates.StsPumpRatesWizard", "Pump Rates",
    		        "Load pump rate data"	,
    		        "Must have opened or created a project.",
    		        "PumpRates20x20.gif");

		    StsWorkflowTreeNode loadTankLevels = loadGroup.addChild("com.Sts.Actions.Wizards.FlowSystem.TankLevels.StsTankLevelsWizard", "Tank Levels",
	    	        "Load storage tank levels",
	    	        "Must have opened or created a project.",
	    	        "TankLevels20x20.gif");

            treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, loadPumpRates));
            treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, loadTankLevels));
            treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, loadDtsRates));
            treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, loadMicroseismicPrestack));
            treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, openSpiritWells));

        StsWorkflowTreeNode defineGroup = workflowRoot.addChild("Define / Edit", null);
            StsWorkflowTreeNode flowNodes = defineGroup.addChild("com.Sts.Actions.Wizards.FlowSystem.FlowNodes.StsFlowNodesWizard", "Flow Nodes",
                    "Interactively track fractures through microseismic events.",
                    "Must have loaded at least one sensor dataset. This step allows the interactive tracking of sensor events.",
                    "FlowNodes20x20.gif");
            nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { loadPumpRates, loadTankLevels}, StsNodeBundle.ALL_REQUIRED);
            treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, flowNodes));

        StsWorkflowTreeNode analysisGroup = workflowRoot.addChild("Analyze", null);
        	StsWorkflowTreeNode fractureTracker = analysisGroup.addChild("com.Sts.Actions.Wizards.FractureTrack.StsFractureTrackWizard", "Fracture Tracking",
                	"Interactively track fractures through microseismic events.",
                	"Must have loaded at least one sensor dataset. This step allows the interactive tracking of sensor events.",
                	"fractureTracker20x20.gif");
        	StsWorkflowTreeNode fractureInterpret = analysisGroup.addChild("com.Sts.Actions.Wizards.FractureInterpret.StsFractureInterpretWizard", "Interpret Fractures",
                	"Interpret clustered events into fractures and fracture sets.",
                	"Must have loaded and clustered at least one sensor dataset. This step allows the interactive tracking of sensor events.",
                	"fractureInterpret20x20.gif");
        	StsWorkflowTreeNode surfaceCurvature = analysisGroup.addChild("com.Sts.Actions.Wizards.SurfaceCurvature.StsSurfaceCurvatureWizard", "Surface Curvature",
                	"Compute Curvature Attributes on a Surface",
                	"Must have loaded at least one surface.",
                	"importSurfaces20x20.gif");
        	StsWorkflowTreeNode volumeCurvature = analysisGroup.addChild("com.Sts.Actions.Wizards.SurfaceCurvature.StsVolumeCurvatureWizard", "Volume Curvature",
                	"Compute Curvature Attributes on a Volume",
                	"Must have loaded at least one seismic volume.",
                	"loadSeismic20x20.gif");
        	StsWorkflowTreeNode logToVolume = analysisGroup.addChild("com.Sts.Actions.Wizards.LogToVolume.StsLogToVolumeWizard", "Log To Volume",
                	"Compute a Volume from Log Data",
                	"Must have loaded at least one well with logs.",
                	"well20x20.gif");
        	StsWorkflowTreeNode sensorCorrelation = analysisGroup.addChild("com.Sts.Actions.Wizards.MicroseismicCorrelation.StsMicroseismicCorrelationWizard", "Sensor Seismic Correlation",
                	"Correlate Sensors with Seismic Attribute",
                	"Must have loaded at least one seismic volume and one dynamic sensor.",
                	"fractureInterpret20x20.gif");

            treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, fractureInterpret));
            treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, fractureTracker));
            treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, volumeCurvature));
            treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, surfaceCurvature));
            treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, sensorCorrelation));

            /*
            nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { getNode(LOAD_DATA, L_POSTSTACK3D), getNode(PROCESS_SEISMIC, P_POSTSTACK3D) }, StsNodeBundle.ONE_REQUIRED);
            treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, volumeCurvature));
            nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { getNode(LOAD_DATA, L_SURFACE) }, StsNodeBundle.ONE_REQUIRED);
            treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, surfaceCurvature));
            */

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
}