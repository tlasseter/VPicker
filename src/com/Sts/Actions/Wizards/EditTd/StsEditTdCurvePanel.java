package com.Sts.Actions.Wizards.EditTd;

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
import com.Sts.MVC.ViewWell.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;
import com.Sts.Types.*;

import javax.swing.*;
import java.awt.*;
import com.Sts.MVC.StsActionPanel;

public class StsEditTdCurvePanel extends StsActionPanel
{
    StsEditTdCurve tdEdit;
    StsWizard wizard;
    //StsFieldBeanPanel wellPanel = new StsFieldBeanPanel();
	StsJPanel wellPanel = StsJPanel.addInsets();
	StsJPanel wellWindowPanel = null;
	StsWellViewModel wellViewModel = null;
	StsWellViewModel prevWellViewModel = null;
    StsComboBoxFieldBean wellListBean = new StsComboBoxFieldBean();
    public StsComboBoxFieldBean tdEditListBean = new StsComboBoxFieldBean();

    StsWell well = null;

//    public static int displayableHeight = 600;

    public StsEditTdCurvePanel(StsEditTdCurve editTDCurve)
    {
        this.tdEdit = editTDCurve;
        this.wizard = tdEdit.getWizard();
		initializeLists();
        initializeWellPanel();
		//wellPanel.setPreferredSize(StsWellViewModel.getDefaultWellWindowSize());
    }

	public void initialize(Object[] wells)
	 {
		 wellListBean.setListItems(wells);
		 tdEditListBean.setListItems(new Object[] {StsEditTd.getNullTdEdit()});
	 }

	private void initializeLists()
	 {
		 wellListBean.initialize(this, "selectedWell", "Select Well:");
		 tdEditListBean.initialize(tdEdit, "tdEdit", "Current td edit:");
	 }

	public void setSelectedWellInList(Object wellObject)
	{
		wellListBean.setSelectedItem(wellObject);
	}

    public void tdEditAdded(StsEditTd addedTdEdit)
    {
        //wellListBean.initialize(this, "selectedWell", "Selected well:");
    }

    void initializeWellPanel()
    {

	    wellPanel = StsJPanel.addInsets();
		wellPanel.setLayout(new BoxLayout(wellPanel,BoxLayout.Y_AXIS));
		//gbc.fill = gbc.HORIZONTAL;
		gbc.fill = GridBagConstraints.BOTH;

	    add(wellPanel);

	    wellPanel.gbc.fill = gbc.HORIZONTAL;
	    wellPanel.addEndRow(wellListBean);

        wellPanel.addEndRow(tdEditListBean);
        StsButton createButton = new StsButton("Create td edit",
                                               "Create a new td edit which is a copy of the current one.", tdEdit,
                                               "createTdEdit");
        wellPanel.addEndRow(createButton);
        wellPanel.setBorder(BorderFactory.createEtchedBorder());

    }


    public void reshape(int x, int y, int width, int height)
    {
        super.reshape(x, y, width, height);
    }

    public Object getSelectedWell()
    {
        return tdEdit.getSelectedWell();
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

    public void setSelectedWell(Object wellObject)
    {
        well = null;

        try
        {
            setSelectedWellInList(wellObject);
            tdEdit.setSelectedWell(wellObject);
            well = (StsWell) wellObject;
//System.out.println("well is "+well.getName());
            if (wellViewModel != null)
            {
                wellPanel.remove(wellViewModel.wellWindowPanel);
                wellViewModel = null;
            }
			if (!well.getName().equals("none"))
			{
				wellViewModel = new StsWellWizardViewModel(well, wizard);
				System.out.println("well view model = "+wellViewModel+" "+well.wellViewModel);
				actionManager = wellViewModel.actionManager;
				// setToPreviousWellViewModel();
				wellPanel.gbc.fill = GridBagConstraints.VERTICAL;
				//wellPanel.gbc.gridwidth = 2;   // Combobox requires 2 grids width
				wellPanel.addEndRow(wellViewModel.wellWindowPanel);

				wellViewModel.initAndStartGL();

				Object[] tdEditsList = tdEdit.getEditTdSet(well).getTdEdits();
				tdEditListBean.setListItems(tdEditsList);
				tdEditListBean.setToLastItem();
				rebuild();
			}
        }
        catch (Exception e)
        {
			e.printStackTrace();
            StsException.outputException("StsEditTdCurvePanel.setSelectedWell() failed for well " + well.getName(),
                                         e, StsException.WARNING);
        }
    }

    public void initializeTdEditListBean()
    {
        Object[] tdEditsList = tdEdit.getEditTdSet(well).getTdEdits();
        tdEditListBean.removeAll();
        tdEditListBean.setListItems(tdEditsList);
        tdEditListBean.setToLastItem();
    }

    public void setTdEdit(Object tdEdit)
	{
		this.tdEdit = (StsEditTdCurve) tdEdit;

    }

	public Object getTdEdit()
	{
	    return tdEdit;
    }
    public void rebuild()
    {
        wizard.rebuild();
    }

}
