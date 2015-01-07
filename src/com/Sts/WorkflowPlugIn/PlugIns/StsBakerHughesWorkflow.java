package com.Sts.WorkflowPlugIn.PlugIns;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.Utilities.*;
import com.Sts.Workflow.*;
import com.Sts.WorkflowPlugIn.*;

public class StsBakerHughesWorkflow extends StsWorkflowPlugIn
{
    public StsBakerHughesWorkflow()
    {
        name = "StsBakerHughesWorkflow";
        workflowName = "Baker Hughes Workflow";
        checkName();
        description = new String("This workflow is used to optimize pump rates to " +
        		"maximize tank capacity. Pump rates and tank levels are monitored and analyzed, " +
        		" rates are adjusted to maximize tank usage through optimal pump rates.");
        createComboBoxDescriptors();
    }

    public void createComboBoxDescriptors()
    {
        ComboBoxDescriptor comboBoxDescriptor;
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

        createGroupAndNodes(workflowRoot, treeModel, LOAD_DATA, new byte[] {L_WELLSANDLOGS, L_SENSOR});
		StsWorkflowTreeNode loadDtsRates = getGroupNode(LOAD_DATA).addChild("com.Sts.Actions.Wizards.DTS.BhDtsLoadWizard", "DTS Files",
    		"Load DtsFiles."	,
    		"Must have loaded wells to associate with this DTS dataset.",
    		"pumpRates20x20.gif");

		StsWorkflowTreeNode loadPumpRates = getGroupNode(LOAD_DATA).addChild("com.Sts.Actions.Wizards.FlowSystem.PumpRates.StsPumpRatesWizard", "Pump Rates",
    		"Load pump rate data for analysis."	,
    		"Must have opened or created a project.",
    		"pumpRates20x20.gif");

		StsWorkflowTreeNode loadTankLevels = getGroupNode(LOAD_DATA).addChild("com.Sts.Actions.Wizards.FlowSystem.TankLevels.StsTankLevelsWizard", "Tank Levels",
	    	"Load tank levels for analysis.",
	    	"Must have opened or created a project.",
	    	"tankLevels20x20.gif");

        StsWorkflowTreeNode flowNodes = getGroupNode(LOAD_DATA).addChild("com.Sts.Actions.Wizards.FlowSystem.FlowNodes.StsFlowNodesWizard", "Flow Nodes",
            "Interactively track fractures through microseismic events.",
            "Must have loaded at least one sensor dataset. This step allows the interactive tracking of sensor events.",
            "flowNodes20x20.gif");

        nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { getNode(LOAD_DATA, L_WELLSANDLOGS)}, StsNodeBundle.ONE_REQUIRED);
        treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, loadDtsRates));
        treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, loadPumpRates));
        treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, loadTankLevels));
        nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { getNode(LOAD_DATA, L_SENSOR) }, StsNodeBundle.ONE_REQUIRED);
        treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, flowNodes));
        
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
            StsTreeNode rootNode = objectTree.createRootNode(model.getProject(), "Project");

            StsTreeNode dataNode = rootNode.addStaticNode("Data");

            dataNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsTank"), "Storage Tanks");
            dataNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsPump"), "Pumps");
            StsSensorClass sensorClass = (StsSensorClass)model.getCreateStsClass("com.Sts.DBTypes.StsSensor");

            StsTreeNode sensorNode = dataNode.addDynamicNode(sensorClass, "Sensor");
			StsSubType[] subTypes = sensorClass.getSubTypes();
			for(int n = 0; n < subTypes.length; n++)
				sensorNode.addDynamicNode(subTypes[n], subTypes[n].getName());

            dataNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsWell"), "Wells");
            dataNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsLogCurveType"), "Log Types");
            dataNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsTimeLogCurveType"), "Time Log Types");
            dataNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsSpectrum"), "Palettes");
            dataNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsCultureObjectSet2D"), "Culture Sets");
            dataNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsMovie"), "Movies");
            dataNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsMonitor"), "Monitored");
            dataNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsWeighPoint"), "WayPoints");
            objectTree.finalizeTreeModel();
            return true;

        }
        catch(Exception e)
        {
            StsException.outputException("StsWorkflowPlug.createObjectsPanel() failed.",
                e, StsException.WARNING);
            return false;
        }
    }
}
