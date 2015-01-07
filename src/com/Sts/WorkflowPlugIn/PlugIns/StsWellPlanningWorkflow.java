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

public class StsWellPlanningWorkflow extends StsWorkflowPlugIn
{
    public StsWellPlanningWorkflow()
    {
        name = "StsWellPlanningWorkflow";
        workflowName = "Platform and Well Planning";
        checkName();
        description = new String("The Well Planning Workflow provides the tools to define platforms" +
                                 " with slot configurations and then assign existing and planned" +
                                 " wells to the slots. Additionally, the workflow includes well planning." +
                                 " A velocity model must exist prior to planning a well.");
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

        //StsWorkflowTreeNode workflowRoot = root.addChild("Platform & Well Planning", null);
        
        createGroupAndNodes(workflowRoot, treeModel, PROCESS_SEISMIC, new byte[] {P_POSTSTACK2D, P_POSTSTACK3D, P_VSP});        
        createGroupAndNodes(workflowRoot, treeModel, LOAD_DATA, new byte[] {L_POSTSTACK2D, L_POSTSTACK3D, L_WELLSANDLOGS, L_VSP, L_COLORPALETTE, L_CULTURE, L_ANCILLARY, L_SURFACE});        

        // Custom addition to Load Group
        StsWorkflowTreeNode workflowLoadWellPlan = getGroupNode(LOAD_DATA).addChild("com.Sts.Actions.Wizards.WellPlan.StsLoadWellPlanWizard", "Well Plan",
            "Load a well plan in UT or S2S format into the project.",
            "A velocity model must exist prior to loading a well plan. Proceed to the Velocity Model Workflow to construct velocity model.",        		
        	"drillPlan20x20.gif");                   
        
        StsWorkflowTreeNode defineData = workflowRoot.addChild("Define / Edit", null);        	
        StsWorkflowTreeNode workflowDefinePlatform = defineData.addChild("com.Sts.Actions.Wizards.PlatformPlan.StsPlatformPlanWizard", "Platform",
            "Define a platform with type, angle, slot configuration and location. Platform will be used to define future well plans.",
            "A project must be open or created prior to defining a platform.",
            "drillPlan20x20.gif");
        StsWorkflowTreeNode workflowDefineWellPlan = defineData.addChild("com.Sts.Actions.Wizards.WellPlan.StsDefineWellPlanWizard", "Well Plan",
            "Define a well plan with type, location, kick-off points and seismic targets that conforms to engineering limits. Well plans will be automatically revised in depth as velocity models are adjusted.",
            "A project must be open or created prior to defining a platform.",
            "drillPlan20x20.gif");
        StsWorkflowTreeNode workflowDefineHorizons = defineData.addChild("com.Sts.Actions.Wizards.Horizons.StsHorizonsWizard", "Horizons",
            "Horizon definition is used to ensure that all surfaces are translated to the same grid and optionally related to markers. It is required for fault definition and reservoir model construction.",
            "Must have loaded or auto-picked at least one surface. Either run the Load->Surfaces or Define->Pick Surfaces workflow step(s) prior to defining horizons.",
            "defineHorizons20x20.gif");
        StsWorkflowTreeNode workflowVelocityModel = defineData.addChild("com.Sts.Actions.Wizards.Velocity.StsVelocityWizard", "Velocity Model",
            "Construct a post-stack velocity model using any combination of user input, and imported well markers, horizons and velocity data.",
            "The minimum reqiurement to construct a velocity model is horizon data or an imported velocity cube. Run the Define->Horizons and/or Load->3D PostStack Volume first. If well ties are desired, import well markers with wells via the Load->Wells & Logs workflow step.",
            "loadSeismic20x20.gif");
            
            //nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { workflowVelocityModel }, StsNodeBundle.ONE_REQUIRED);
            //treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, workflowDefineWellPlan));
            treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, workflowDefinePlatform));
            treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, workflowLoadWellPlan));
            treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, workflowDefineWellPlan));
            nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { workflowDefineHorizons, getNode(LOAD_DATA, L_POSTSTACK3D), getNode(LOAD_DATA, L_SURFACE)}, StsNodeBundle.ONE_REQUIRED);
            treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, workflowDefineHorizons));
            treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, workflowVelocityModel));
            //nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { workflowVelocityModel }, StsNodeBundle.ONE_REQUIRED);            
            treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, workflowLoadWellPlan));
/*            
        StsWorkflowTreeNode miscSteps = workflowRoot.addChild("Output", null);
          StsWorkflowTreeNode makeMovie = miscSteps.addChild("com.Sts.Actions.Wizards.MakeMovie.StsMakeMovieWizard", "QT Movie", "loadAncillaryData20x20.gif");
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
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsSeismicLineSet"), "2D Line Sets", true);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsPreStackLineSet3d"), "3D Gathers", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsPreStackLineSet2d"), "2D Gathers", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsPreStackMicroseismicSet"), "Microseismic Gathers", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsLiveWell"), "Live Wells", false);                        
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsWell"), "Wells", true);
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
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsPlatform"), "Drilling Platform", true);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsWellPlanSet"), "Well Plans", true);

            StsTreeNode sensorNode = checkAddStaticNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsSensor"), "Sensors", false);
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
            StsException.outputException("StsFractureAnalysisWorkflow.createObjectsPanel() failed.",
                e, StsException.WARNING);
            return false;
        }
}
}
