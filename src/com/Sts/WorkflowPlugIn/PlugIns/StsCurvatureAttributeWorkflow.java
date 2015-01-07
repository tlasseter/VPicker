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

public class StsCurvatureAttributeWorkflow extends StsWorkflowPlugIn
{
    public StsCurvatureAttributeWorkflow()
    {
        name = "StsCurvatureAttributeWorkflow";
        workflowName = "Seismic Attributes";
        checkName();
        description = new String("Calculates curvature attribute on surfaces and seismic volumes." );

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

        //createGroupAndNodes(workflowRoot, treeModel, PROCESS_SEISMIC, new byte[] {P_PRESTACK2D, P_PRESTACK3D, P_POSTSTACK2D, P_POSTSTACK3D, P_VSP});
        //createGroupAndNodes(workflowRoot, treeModel, LOAD_DATA, new byte[] {L_PRESTACK2D, L_PRESTACK3D, L_POSTSTACK2D, L_POSTSTACK3D, L_WELLSANDLOGS, L_VSP, L_SENSOR, L_SURFACE, L_COLORPALETTE, L_CULTURE, L_ANCILLARY});
        createGroupAndNodes(workflowRoot, treeModel, LOAD_DATA, new byte[] {L_POSTSTACK3D, L_SURFACE});

        StsWorkflowTreeNode defineData = workflowRoot.addChild("Define / Edit", null);
//        	StsWorkflowTreeNode virtualVolume = defineData.addChild("com.Sts.Actions.Wizards.VirtualVolume.StsVirtualVolumeWizard", "Virtual Volumes",
//        		"Create a virtual volume. Virtual volumes are seismic volumes that are computed on the fly from other seismic attribute volumes or virtual volumes. Filter volumes, crossplot, blended and math volumes are all virtual volume types. Virtual volumes are treated as real volumes everywhere else in S2S."	,
//        		"Must have at least one seismic volume to define a virtual volume. run the Load->3D PostStack Seismic workflow step first.",
//        		"defineVirtual20x20.gif");
        	StsWorkflowTreeNode volumeFilter = defineData.addChild("com.Sts.Actions.Wizards.VolumeFilter.StsVolumeFilterWizard", "Volume Filter",
        		"Run enhanced filters on existing seismic volumes.",
        		"Must have at least one seismic volume loaded to run the seismic attribute step. Run the Load->3D Poststack Seismic Workflow step first.",
        		"seismicAttribute20x20.gif");
//        	//StsWorkflowTreeNode interpolateVolume = defineData.addChild("com.Sts.Actions.Wizards.VolumeInterpolation.StsVolumeInterpolationWizard", "Interpolate Volume",
//            	//"Interpolate null values and missing traces from attribute volumes.",
//            	//"Must have at least one seismic volume loaded to run the interpolate volume step. Run the Load->3D Poststack Seismic Workflow step first.",
//            	//"seismicAttribute20x20.gif");
//        	StsWorkflowTreeNode subVolume = defineData.addChild("com.Sts.Actions.Wizards.SubVolume.StsSubVolumeWizard", "SubVolumes",
//                "Define a sub-volume of the current project. Sub-volumes can be defined in relation to surfaces or as an arbitrary set of polyhedron. Sub-volumes can be applied to any seismic volume or crossplot to limit the amount of data displayed on the cursor planes or crossplotted."		,
//                "A minimum of one seismic volume must exist prior to defining a subvolume. Run the Load->3D Poststack Seismic before defining a subvolume.",
//        		"defineCropVolume20x20.gif");
//        	StsWorkflowTreeNode subVolumeMgmt = defineData.addChild("com.Sts.Actions.Wizards.SubVolumeMgmt.StsSubVolumeMgmtWizard", "Manage SubVolumes",
//        		"Manage the application of sub-volumes to the various data elements. This is simply a convenience function and everything that can be done on this workflow step can also be done from the object panel.",
//                "Must have defined at least one sub-volume prior to attempting to manage sub-volumes. Seems obvious. Run the Define->SubVolume workflow step first.",
//        		"defineCropVolume20x20.gif");
//        	StsWorkflowTreeNode movie = defineData.addChild("com.Sts.Actions.Wizards.Movie.StsMovieWizard", "Animations",
//            	"Animate the cursor planes and user views. control the speed and interval a slice is moved through the volume while rotating the users azimuth and elevation relative to the data.",
//            	"No pre-requisites to create and run an animation. May not be very interesting without data, but it will work.",
//        		"createMovie20x20.gif");

        	StsNodeBundle nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { getNode(LOAD_DATA, L_POSTSTACK3D), getNode(PROCESS_SEISMIC, P_POSTSTACK3D) }, StsNodeBundle.ONE_REQUIRED);
        	treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, volumeFilter));
//        	treeModel.addNodeConnection(new StsNodeConnection(projectNodeBundle, movie));
//            treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, seismicAttribute));
//            //treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, interpolateVolume));
//            treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, subVolume));
//            treeModel.addNodeConnection(new StsNodeConnection(subVolume, subVolumeMgmt));
            
            StsWorkflowTreeNode analysisGroup = workflowRoot.addChild("Analyze", null);
        	StsWorkflowTreeNode surfaceCurvature = analysisGroup.addChild("com.Sts.Actions.Wizards.SurfaceCurvature.StsSurfaceCurvatureWizard", "Surface Curvature",
                	"Compute Curvature Attributes on a Surface",
                	"Must have loaded at least one surface.",
                	"importSurfaces20x20.gif");
        	StsWorkflowTreeNode volumeCurvature = analysisGroup.addChild("com.Sts.Actions.Wizards.SurfaceCurvature.StsVolumeCurvatureWizard", "Volume Curvature",
                	"Compute Curvature Attributes on a Volume",
                	"Must have loaded at least one seismic volume.",
                	"loadSeismic20x20.gif");
//        	nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { getNode(LOAD_DATA, L_SURFACE)}, StsNodeBundle.ONE_REQUIRED);
//        	treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, surfaceCurvature));
        	nodeBundle = new StsNodeBundle(new StsWorkflowTreeNode[] { getNode(LOAD_DATA, L_POSTSTACK3D) }, StsNodeBundle.ONE_REQUIRED);
        	treeModel.addNodeConnection(new StsNodeConnection(nodeBundle, volumeCurvature));
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
}
