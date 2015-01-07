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

public class StsOpenSpiritWorkflow extends StsWorkflowPlugIn
{
    public StsOpenSpiritWorkflow()
    {
        name = "StsOpenSpiritWorkflow";
        workflowName = "OpenSpirit";
        checkName();
        description = new String("Connects to OpenSpirit Server and loads data through the connection from all" +
        		" OpenSpirit supportted datastores.");

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

        StsWorkflowTreeNode defineData = workflowRoot.addChild("Connect", null);
        	StsWorkflowTreeNode openSpiritConnect = defineData.addChild("com.Sts.Actions.Wizards.OSConnect.StsOSConnectWizard", "OpenSpirit",
        		"Connect to an OpenSpirit Data Server.",
        		"A project must be open prior to running this workflow step.",
        		"connectOpenSpirit20x20.gif");
        	
        StsWorkflowTreeNode loadGroup = workflowRoot.addChild("Load OpenSpirit", null);
    		StsWorkflowTreeNode openSpiritWells = loadGroup.addChild("com.Sts.Actions.Wizards.OSWell.StsOSWellWizard", "Wells",
        		"Load well and well related data from an OpenSpirit datastore."	,
        		"Must have established an OpenSpirit connection prior to attempting to load data.",
        		"openSpiritWells20x20.gif"); 
        	StsNodeBundle nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { openSpiritConnect }, StsNodeBundle.ONE_REQUIRED);        	
            treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, openSpiritWells));    		
        
       logUsageChange();
    }

    public void addAdditionalToolbars(StsWin3dBase win3d)
    {
    	;
    }

    public void addOptionalNodes(StsTreeModel treeModel, String[] options)
    {
        ;
    }
}
