package com.Sts.WorkflowPlugIn.PlugIns;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.Utilities.*;
import com.Sts.Workflow.*;
import com.Sts.WorkflowPlugIn.*;

public class StsPumpTankWorkflow extends StsWorkflowPlugIn
{
    public StsPumpTankWorkflow()
    {
        name = "StsPumpTankWorkflow";
        workflowName = "Pump/Tank Optimization Workflow";
        checkName();
        description = new String("This workflow is used to optimize pump rates to " +
        		"maximize tank capacity. Pump rates and tank levels are monitored and analyzed, " +
        		" rates are adjusted to maximize tank usage through optimal pump rates.");
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

        createGroupAndNodes(workflowRoot, treeModel, LOAD_DATA, new byte[] {L_WELLSANDLOGS});               
		StsWorkflowTreeNode loadPumpRates = getGroupNode(LOAD_DATA).addChild("com.Sts.Actions.Wizards.FlowSystem.PumpRates.StsPumpRatesWizard", "Load Pump Rates",
    		"Load pump rate data for analysis."	,
    		"Must have opened or created a project.",
    		"pumpRates20x20.gif");
		
		StsWorkflowTreeNode loadTankLevels = getGroupNode(LOAD_DATA).addChild("com.Sts.Actions.Wizards.FlowSystem.TankLevels.StsTankLevelsWizard", "Load Tank Levels",
	    		"Load tank levels for analysis.",
	    		"Must have opened or created a project.",
	    		"tankLevels20x20.gif");

        	StsWorkflowTreeNode flowNodes = getGroupNode(LOAD_DATA).addChild("com.Sts.Actions.Wizards.FlowSystem.FlowNodes.StsFlowNodesWizard", "Flow Nodes",
                	"Interactively track fractures through microseismic events.",
                	"Must have loaded at least one sensor dataset. This step allows the interactive tracking of sensor events.",
                	"flowNodes20x20.gif"); 

        	nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { loadPumpRates, loadTankLevels }, StsNodeBundle.ALL_REQUIRED);        	
            treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, loadPumpRates));
            treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, loadTankLevels)); 
//            treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, chokePoints));            
        
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
            dataNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsPlatform"), "Drilling Platform");
            dataNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsWellPlanSet"), "Well Plans");
            dataNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsWell"), "Wells");
            dataNode.addDynamicNode(model.getCreateStsClass("com.Sts.DBTypes.StsLogCurveType"), "Log Types");

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
