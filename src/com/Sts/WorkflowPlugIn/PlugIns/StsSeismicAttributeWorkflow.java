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

public class StsSeismicAttributeWorkflow extends StsWorkflowPlugIn
{
    public StsSeismicAttributeWorkflow()
    {
        name = "StsSeismicAttributeWorkflow";
        workflowName = "Seismic Attributes";
        checkName();
        description = new String("Loads, filters and calculates seismic volumes, and displays volumes," +
        		" with well, logs, sensors and surfaces");

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

        createGroupAndNodes(workflowRoot, treeModel, PROCESS_SEISMIC, new byte[] {P_PRESTACK2D, P_PRESTACK3D, P_POSTSTACK2D, P_POSTSTACK3D, P_VSP});
        createGroupAndNodes(workflowRoot, treeModel, LOAD_DATA, new byte[] {L_PRESTACK2D, L_PRESTACK3D, L_POSTSTACK2D, L_POSTSTACK3D, L_WELLSANDLOGS, L_VSP, L_SENSOR, L_SURFACE, L_COLORPALETTE, L_CULTURE, L_ANCILLARY});

        StsWorkflowTreeNode defineData = workflowRoot.addChild("Define / Edit", null);
        	StsWorkflowTreeNode virtualVolume = defineData.addChild("com.Sts.Actions.Wizards.VirtualVolume.StsVirtualVolumeWizard", "Virtual Volumes",
        		"Create a virtual volume. Virtual volumes are seismic volumes that are computed on the fly from other seismic attribute volumes or virtual volumes. Filter volumes, crossplot, blended and math volumes are all virtual volume types. Virtual volumes are treated as real volumes everywhere else in S2S."	,
        		"Must have at least one seismic volume to define a virtual volume. run the Load->3D PostStack Seismic workflow step first.",
        		"defineVirtual20x20.gif");
        	StsWorkflowTreeNode seismicAttribute = defineData.addChild("com.Sts.Actions.Wizards.SeismicAttribute.StsSeismicAttributeWizard", "Seismic Attributes",
        		"Compute seismic attributes from existing seismic volumes.",
        		"Must have at least one seismic volume loaded to run the seismic attribute step. Run the Load->3D Poststack Seismic Workflow step first.",
        		"seismicAttribute20x20.gif");
        	//StsWorkflowTreeNode interpolateVolume = defineData.addChild("com.Sts.Actions.Wizards.VolumeInterpolation.StsVolumeInterpolationWizard", "Interpolate Volume",
            	//"Interpolate null values and missing traces from attribute volumes.",
            	//"Must have at least one seismic volume loaded to run the interpolate volume step. Run the Load->3D Poststack Seismic Workflow step first.",
            	//"seismicAttribute20x20.gif");
        	StsWorkflowTreeNode subVolume = defineData.addChild("com.Sts.Actions.Wizards.SubVolume.StsSubVolumeWizard", "SubVolumes",
                "Define a sub-volume of the current project. Sub-volumes can be defined in relation to surfaces or as an arbitrary set of polyhedron. Sub-volumes can be applied to any seismic volume or crossplot to limit the amount of data displayed on the cursor planes or crossplotted."		,
                "A minimum of one seismic volume must exist prior to defining a subvolume. Run the Load->3D Poststack Seismic before defining a subvolume.",
        		"defineCropVolume20x20.gif");
        	StsWorkflowTreeNode subVolumeMgmt = defineData.addChild("com.Sts.Actions.Wizards.SubVolumeMgmt.StsSubVolumeMgmtWizard", "Manage SubVolumes",
        		"Manage the application of sub-volumes to the various data elements. This is simply a convenience function and everything that can be done on this workflow step can also be done from the object panel.",
                "Must have defined at least one sub-volume prior to attempting to manage sub-volumes. Seems obvious. Run the Define->SubVolume workflow step first.",
        		"defineCropVolume20x20.gif");
        	StsWorkflowTreeNode movie = defineData.addChild("com.Sts.Actions.Wizards.Movie.StsMovieWizard", "Animations",
            	"Animate the cursor planes and user views. control the speed and interval a slice is moved through the volume while rotating the users azimuth and elevation relative to the data.",
            	"No pre-requisites to create and run an animation. May not be very interesting without data, but it will work.",
        		"createMovie20x20.gif");

        	StsNodeBundle nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { getNode(LOAD_DATA, L_POSTSTACK3D), getNode(PROCESS_SEISMIC, P_POSTSTACK3D) }, StsNodeBundle.ONE_REQUIRED);
        	treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, virtualVolume));
        	treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, movie));
            treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, seismicAttribute));
            //treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, interpolateVolume));
            treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, subVolume));
            treeModel.addNodeConnection(new StsNodeConnection(subVolume, subVolumeMgmt));
            
            StsWorkflowTreeNode analysisGroup = workflowRoot.addChild("Analyze", null);
        	StsWorkflowTreeNode multiAttributeAnalysis = analysisGroup.addChild("com.Sts.Actions.Wizards.MultiAttributeAnalysis.StsMultiAttrAnalysisWizard", "3D Attribute Vector",
        		"Create a dynamic map of vectors on the Z slice or on a surface, where the azimuth, length and color are controlled by three seismic volumes. This capability is often used to analyze natural fracture patterns.",
        		"Must have loaded or created at least two seismic volumes to run the 3D Vector Workflow Step. Run Load->3D PostStack Seismic or Define->Seismic Attributes (Seismic Attribute Workflow) prior to running 3D Attribute Vector",
        		"loadCulture20x20.gif");
        	nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { getNode(LOAD_DATA, L_POSTSTACK3D), getNode(PROCESS_SEISMIC, P_POSTSTACK3D) }, StsNodeBundle.ONE_REQUIRED);
        	treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, multiAttributeAnalysis));
/*
        StsWorkflowTreeNode outputData = workflowRoot.addChild("Output", null);
          StsWorkflowTreeNode makeMovie = outputData.addChild("com.Sts.Actions.Wizards.MakeMovie.StsMakeMovieWizard", "QT Movie", "loadAncillaryData20x20.gif");
          treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, makeMovie));
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
        ;
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
            StsTreeNode virtualVolumeNode = checkAddStaticNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsVirtualVolume"), "Virtual Volumes", true);
            if(virtualVolumeNode != null)
            {
                checkAddDynamicNode(model, virtualVolumeNode, model.getCreateStsClass("com.Sts.DBTypes.StsMathVirtualVolume"), "Math", true);
                checkAddDynamicNode(model, virtualVolumeNode, model.getCreateStsClass("com.Sts.DBTypes.StsBlendedVirtualVolume"), "Blended", true);
                checkAddDynamicNode(model, virtualVolumeNode, model.getCreateStsClass("com.Sts.DBTypes.StsCrossplotVirtualVolume"), "Crossplot", true);
                checkAddDynamicNode(model, virtualVolumeNode, model.getCreateStsClass("com.Sts.DBTypes.StsRGBAVirtualVolume"), "RGBA", false);
                checkAddDynamicNode(model, virtualVolumeNode, model.getCreateStsClass("com.Sts.DBTypes.StsFilterVirtualVolume"), "Filter", true);
                checkAddDynamicNode(model, virtualVolumeNode, model.getCreateStsClass("com.Sts.DBTypes.StsSensorVirtualVolume"), "Sensor", true);
            }
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsCrossplot"), "Crossplot", false);

            subVolumeNode = checkAddStaticNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsSubVolume"), "SubVolumes", true);
            if(subVolumeNode != null)
            {
                checkAddDynamicNode(model, subVolumeNode, model.getCreateStsClass("com.Sts.DBTypes.StsDualSurfaceSubVolume"), "Dual Surface", true);
                checkAddDynamicNode(model, subVolumeNode, model.getCreateStsClass("com.Sts.DBTypes.StsBoxSetSubVolume"), "Box Set", true);
                checkAddDynamicNode(model, subVolumeNode, model.getCreateStsClass("com.Sts.DBTypes.StsWellSubVolume"), "Well Set", false);
            }
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsHorpick"), "Horizon Picks", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsSurface"), "Surfaces", true);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsPlatform"), "Drilling Platform", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsWellPlanSet"), "Well Plans", false);

            StsTreeNode sensorNode = checkAddStaticNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsSensor"), "Sensors", true);
            checkAddDynamicNode(model, sensorNode, model.getCreateStsClass("com.Sts.DBTypes.StsStaticSensor"), "Static", true);
            checkAddDynamicNode(model, sensorNode, model.getCreateStsClass("com.Sts.DBTypes.StsDynamicSensor"), "Dynamic", true);

            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsMVFractureSet"), "MV Fractures", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsTriangulatedFracture"), "Fractures", false);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsMultiAttributeVector"), "Multi-Attribute Vectors", true);
            checkAddDynamicNode(model, dataNode, model.getCreateStsClass("com.Sts.DBTypes.StsFaultStickSet"), "Fault Sticks", false);
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
