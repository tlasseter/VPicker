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

public class StsPreStackVelocityWorkflow extends StsWorkflowPlugIn
{
	private StsTreeModel treeModel;
	public StsPreStackVelocityWorkflow()
	{
		name = "StsPreStackVelocityWorkflow";
		workflowName = "Pre-Stack Velocity Modeling";
		checkName();
		description = new String("The Pre-Stack Velocity Modeling Workflow contains all the steps required" +
				" to analyze and model velocities from migrated and unmigrated" +
		" seismic gathers.\n");
		createComboBoxDescriptors();
	}

	public void createComboBoxDescriptors()
	{
		StsWorkflowPlugIn.ComboBoxDescriptor comboBoxDescriptor;
		String[] volumeClasses = new String[]
		                                    {"com.Sts.DBTypes.StsSeismicVolume",
		                                    "com.Sts.DBTypes.StsFilterVirtualVolume" };
		comboBoxDescriptor =
			constructComboBoxToolbarDescriptor("com.Sts.DBTypes.StsSeismicVolume",
					volumeClasses, "poststack", "noseismic");
		comboBoxToolbarDescriptors.add(comboBoxDescriptor);
		String[] preStackClasses = new String[]
		                                      {"com.Sts.DBTypes.StsPreStackLineSet3d",
		                                      "com.Sts.DBTypes.StsPreStackLineSet2d"};
		comboBoxDescriptor =
			StsWorkflowPlugIn.constructComboBoxToolbarDescriptor("com.Sts.DBTypes.StsPreStackLineSet",
					preStackClasses, "prestack", "noseismic");
		comboBoxToolbarDescriptors.add(comboBoxDescriptor);
	}

	public void createWorkflowNodes(StsTreeModel treeModel,
			StsWorkflowTreeNode workflowRoot)
	{
		StsNodeBundle nodeBundle;
		this.treeModel = treeModel;

		treeModel.addMenuNode(newProject);
		treeModel.addMenuNode(openProject);

		createGroupAndNodes(workflowRoot, treeModel, PROCESS_SEISMIC, new byte[] {P_PRESTACK2D, P_PRESTACK3D, P_POSTSTACK2D, P_POSTSTACK3D, P_VSP});
		createGroupAndNodes(workflowRoot, treeModel, LOAD_DATA, new byte[] {L_PRESTACK2D, L_PRESTACK3D, L_POSTSTACK2D, L_POSTSTACK3D, L_HANDVEL, L_COLORPALETTE});

		StsWorkflowTreeNode defineData = workflowRoot.addChild("Define / Edit",
				null);
		StsWorkflowTreeNode velocityAnalysis =
			defineData.addChild("com.Sts.Actions.Wizards.VelocityAnalysis.StsVelocityAnalysisWizard",
					"Velocity Model", "loadSeismic20x20.gif");
		StsWorkflowTreeNode superGather =
			defineData.addChild("com.Sts.Actions.Wizards.SuperGather.StsSuperGatherWizard",
					"Super Gather",
					"Define a super gather which would then be used for all semblance and CVS displays",
					"Must have a pre-stack velocity model defined prior to defining a super gather. Run the Define->Velocity Model workflow step prior to defining the super gather.",
			"loadSeismic20x20.gif");
		StsWorkflowTreeNode copyVelocity =
			defineData.addChild("com.Sts.Actions.Wizards.CopyVelocity.StsCopyVelocityWizard",
					"Copy Velocity Model",
					"Copy an existing pre-stack velocity model, effectively taking a snap-shot of the current model and allowing continued refinement of copy.",
					"Must have a pre-stack velocity model defined prior to copying it. Run the Define->Velocity Model and pick velocities to create a new model prior to copying it.",
			"loadSeismic20x20.gif");
		StsWorkflowTreeNode filterVolume =
			defineData.addChild("com.Sts.Actions.Wizards.VirtualVolume.StsFilterVolumeWizard",
					"Filter Volume",
					"Filter an existing post-stack volume. Filters primarily consist of smoothing operators with user defined kernal types and sizes. Filter volumes are a type of virtual volume.",
					"Must have a post-stack seismic volume in order to run the filter volume workflow step. Poststack volume can be directly loaded or generated by building a veloicty modeling.",
			"defineVirtual20x20.gif");
		StsWorkflowTreeNode exportData = workflowRoot.addChild("Export",
				null);
		StsWorkflowTreeNode exportPreStack =
			exportData.addChild("com.Sts.Actions.Wizards.PreStackExport.StsPreStackExportWizard",
					"Volumes",
					"Export prestack 2D and 3D data including velocity models, gathers, stacks and semblance.",
					"Must have a pre-stack seismic volume loaded into the project in order to export it or related data elements. Run the Load->PreStack 2D or PreStack 3D Step first.",
			"loadSeismic20x20.gif");
		StsWorkflowTreeNode exportHandVel2D =
			exportData.addChild("com.Sts.Actions.Wizards.PreStackExport.StsPreStackExportHandvelWizard2d",
					"2D Handvel",
					"Export prestack 2D velocity model as Handvel cards.",
					"Must have a pre-stack seismic volume loaded into the project in order to export it or related data elements. Run the Load->PreStack 2D or PreStack 3D Step first.",
			"importHandVels20x20.gif");
		StsWorkflowTreeNode exportHandVel3D =
			exportData.addChild("com.Sts.Actions.Wizards.PreStackExport.StsPreStackExportHandvelWizard3d",
					"3D Handvel",
					"Export prestack 3D velocity model as Handvel cards.",
					"Must have a pre-stack seismic volume loaded into the project in order to export it or related data elements. Run the Load->PreStack 2D or PreStack 3D Step first.",
			"importHandVels20x20.gif");
		nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] {
				getNode(LOAD_DATA, L_PRESTACK3D), getNode(LOAD_DATA,
						L_PRESTACK2D)}, StsNodeBundle.ONE_REQUIRED);
		treeModel.addNodeConnection(new StsNodeConnection(nodeBundle,
				velocityAnalysis));
		treeModel.addNodeConnection(new
				StsNodeConnection(velocityAnalysis, superGather));
		treeModel.addNodeConnection(new
				StsNodeConnection(velocityAnalysis, copyVelocity));
		treeModel.addNodeConnection(new StsNodeConnection(nodeBundle,
				exportPreStack));
		treeModel.addNodeConnection(new StsNodeConnection(nodeBundle,
				exportHandVel2D));
		treeModel.addNodeConnection(new StsNodeConnection(nodeBundle,
				exportHandVel3D));
		nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] {
				getNode(LOAD_DATA, L_POSTSTACK3D) }, StsNodeBundle.ONE_REQUIRED);
		treeModel.addNodeConnection(new StsNodeConnection(nodeBundle,
				filterVolume));

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
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsSeismicLineSet"), "2D Line Sets", true);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsPreStackLineSet3d"), "3D Gathers", true);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsPreStackLineSet2d"), "2D Gathers", true);
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
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsHorpick"), "Horizon Picks", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsSurface"), "Surfaces", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsPlatform"), "Drilling Platform", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsWellPlanSet"), "Well Plans", false);

            StsTreeNode sensorNode = checkAddStaticNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsSensor"), "Sensors", false);
            checkAddDynamicNode(model, sensorNode, model.getCreateStsClass("com.Sts.DBTypes.StsStaticSensor"), "Static", false);
            checkAddDynamicNode(model, sensorNode, model.getCreateStsClass("com.Sts.DBTypes.StsDynamicSensor"), "Dynamic", false);

            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsMVFractureSet"), "MV Fractures", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsTriangulatedFracture"), "Fractures", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsMultiAttributeVector"), "Multi-Attribute Vectors", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsFaultStickSet"), "Fault Sticks", false);
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
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsWeighPoint"), "WayPoints", false);


            // Model Node
            modelNode = rootNode.addStaticNode("Model");
            checkAddDynamicNode(model, modelNode, model.getCreateStsClass("com.Sts.DBTypes.StsModelSurface"), "Horizons", false);
            checkAddDynamicNode(model, modelNode, model.getCreateStsClass("com.Sts.DBTypes.StsPreStackVelocityModel3d"), "PreStack3d 3D Velocity Model", true);
            checkAddDynamicNode(model, modelNode, model.getCreateStsClass("com.Sts.DBTypes.StsPreStackVelocityModel2d"), "PreStack3d 2D Velocity Model", true);
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
            StsException.outputException("StsFractureAnalysisWorkflow.createObjectsPanel() failed.",
                e, StsException.WARNING);
            return false;
        }
}
}