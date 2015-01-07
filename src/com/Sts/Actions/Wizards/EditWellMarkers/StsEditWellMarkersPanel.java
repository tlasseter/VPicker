package com.Sts.Actions.Wizards.EditWellMarkers;

/**
 * <p>Title: jS2S development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2004</p>
 * <p>Company: S2S Systems LLC</p>
 * @author Tom Lasseter
 * @version 1.0
 */

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.ViewWell.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;
import com.Sts.Types.*;
import javax.swing.BoxLayout;

public class StsEditWellMarkersPanel extends StsActionPanel
{
    StsEditWellMarkers editWellMarkers;
    StsWizard wizard = null;
    StsJPanel wellPanel = StsJPanel.addInsets();
    //StsJPanel editMarkerOperationsPanel = null;
    StsComboBoxFieldBean wellListBean = new StsComboBoxFieldBean();
	StsJPanel editMarkersPanel = null;
    StsJPanel wellWindowPanel = null;
    StsWellViewModel wellViewModel = null;
    StsWellViewModel prevWellViewModel = null;
    StsComboBoxFieldBean operationsMenuComboBox = null;
//    JMenuBar editMenuBar = null;
//    JLabel editMenuComboBoxLabel = new JLabel("Selected Marker Action:");
    StsComboBoxFieldBean editMenuComboBox = new StsComboBoxFieldBean();

    StsJPanel wellMarkerActionPanel = null;

    String editType = noActionString; // add or edit
    static public String noActionString = "No Action";
    static public String addMarkerString = "Add Marker";
    static public String editMarkerString = "Edit Marker";

    static public String[] editTypes = new String[]
    {
        noActionString, addMarkerString, editMarkerString
    };

//    public static int displayableHeight = 600;

    public StsEditWellMarkersPanel(StsEditWellMarkers editWellMarkers)
    {
        this.editWellMarkers = editWellMarkers;
        try
        {
            wizard = editWellMarkers.getWizard();
            initializeLists();
            initializeWellPanel();
            wellPanel.setPreferredSize(StsWellViewModel.getDefaultWellWindowSize());
        }
        catch (Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void initialize(Object[] wells)
    {
        wellListBean.setListItems(wells);
    }

    private void initializeLists()
    {
        wellListBean.initialize(this, "selectedWell", "Select Well:");
        editMenuComboBox.initialize(this, "editType", "Select Action:", editTypes);
   }

   public void setSelectedWellInList(Object wellObject)
   {
       wellListBean.setSelectedItem(wellObject);
   }

    private void initializeWellPanel() throws Exception
    {
        wellPanel = StsJPanel.addInsets();
		wellPanel.setLayout(new BoxLayout(wellPanel,BoxLayout.Y_AXIS));
        gbc.fill = gbc.HORIZONTAL;
        addToRow(wellPanel);

        wellPanel.gbc.fill = gbc.HORIZONTAL;
        wellPanel.addEndRow(wellListBean);
    }

    private void initializeEditMarkersPanel() throws Exception
    {
        //editMarkersPanel.addToRow(editMenuComboBoxLabel);
    	editMarkersPanel = StsJPanel.addInsets();
        editMarkersPanel.addEndRow(editMenuComboBox);
        addToRow(editMarkersPanel);
   }

    public void reshape(int x, int y, int width, int height)
    {
        super.reshape(x, y, width, height);
    }

    public Object getSelectedWell()
    {
        return editWellMarkers.getSelectedWell();
    }

    public void setSelectedWell(Object wellObject)
    {
        StsWell well = null;
        try
        {
            setSelectedWellInList(wellObject);
            editWellMarkers.setSelectedWell(wellObject);
            well = (StsWell) wellObject;

            if (wellPanel != null)
			{
				remove(wellPanel);
				wellViewModel=null;
			}
            initializeWellPanel();
            savePreviousWellViewModel();
			if (!well.getName().equals("none"))
            {
                wellViewModel = new StsWellWizardViewModel(well, wizard);
                actionManager = wellViewModel.actionManager;
                // setToPreviousWellViewModel();
                wellPanel.gbc.gridwidth = 2;   // Combobox requires 2 grids width
                wellPanel.addEndRow(wellViewModel.wellWindowPanel);

                if(editMarkersPanel == null)
                	initializeEditMarkersPanel();
                executeAction(editType);
                wellViewModel.initAndStartGL();
            }
            else
            {
                if(editMarkersPanel != null) remove(editMarkersPanel);
                editMarkersPanel = null;
				actionManager=null;
            }
            rebuild();
        }
        catch (Exception e)
        {
            StsException.outputException("StsEditTdCurvePanel.setSelectedWell() failed for well " + well.getName(),
                                         e, StsException.WARNING);
        }
    }

    private void savePreviousWellViewModel()
    {
        if(wellViewModel == null) return;
        prevWellViewModel = wellViewModel;
    }

    private void setToPreviousWellViewModel()
    {
        if(prevWellViewModel == null) return;
    }

    public void rebuild()
    {
        wizard.rebuild();
    }

    public void setEditType(String editType)
    {
        if (this.editType == editType) return;
        this.editType = editType;
        executeAction(editType);
    }

    private void executeAction(String editType)
    {
        if(wellMarkerActionPanel != null)
        	editMarkersPanel.remove(wellMarkerActionPanel);

        // end the current action
		if (actionManager != null)
			actionManager.endCurrentAction();
        editMarkersPanel.gbc.gridwidth = 2;
        if (editType == addMarkerString)
        {
            StsAddWellMarker addWellMarkerAction = new StsAddWellMarker(wizard.model, wellViewModel);
            actionManager.checkStartAction(addWellMarkerAction);
            wellMarkerActionPanel = addWellMarkerAction.getAddMarkerPanel();
            editMarkersPanel.addEndRow(wellMarkerActionPanel);
        }
        else if (editType == editMarkerString)
        {
            if(!wellViewModel.getWell().hasMarkers())
            {
                new StsMessage(wizard.getModel().win3d, StsMessage.WARNING, "No markers exist for this well.");
                return;
            }
            StsEditWellMarker editWellMarkerAction = new StsEditWellMarker(wizard.model, wellViewModel);
            actionManager.checkStartAction(editWellMarkerAction);
            wellMarkerActionPanel = editWellMarkerAction.getEditMarkerPanel();
            editMarkersPanel.addEndRow(wellMarkerActionPanel);
        }

        wizard.rebuild();
    }

    public String getEditType()
    {
        return editType;
    }
}
