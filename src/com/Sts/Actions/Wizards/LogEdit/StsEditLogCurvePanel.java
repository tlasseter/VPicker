package com.Sts.Actions.Wizards.LogEdit;

/**
 * <p>Title: jS2S development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2004</p>
 * <p>Company: S2S Systems LLC</p>
 * @author Tom Lasseter
 * @version 1.0
 */

import com.Sts.Actions.Wizards.StsWizard;
import com.Sts.DBTypes.StsLogCurve;
import com.Sts.DBTypes.StsObjectRefList;
import com.Sts.DBTypes.StsWell;
import com.Sts.DBTypes.StsWellViewModel;
import com.Sts.MVC.StsActionPanel;
import com.Sts.MVC.ViewWell.StsLogCurvesView;
import com.Sts.Types.StsWellWizardViewModel;
import com.Sts.UI.Beans.StsComboBoxFieldBean;
import com.Sts.UI.Beans.StsJPanel;
import com.Sts.UI.StsButton;
import com.Sts.UI.StsToggleButton;
import com.Sts.Utilities.StsException;

import javax.swing.*;
import java.awt.*;

public class StsEditLogCurvePanel extends StsActionPanel
{
    StsLogCurve logCurve;
    StsEditLogCurve editLogCurve;
    StsWizard wizard;
    //StsFieldBeanPanel wellPanel = new StsFieldBeanPanel();
	StsJPanel wellPanel = StsJPanel.addInsets();
	StsJPanel wellWindowPanel = null;
	StsWellViewModel wellViewModel = null;
	StsWellViewModel prevWellViewModel = null;
    StsComboBoxFieldBean wellListBean = new StsComboBoxFieldBean(this, "selectedWell", "Select Well:", new Object[] { nullWell });
    StsComboBoxFieldBean logCurveListBean = new StsComboBoxFieldBean(this, "logCurve", "Current log Curve:", new Object[] { nullLogCurve });
    StsComboBoxFieldBean editBean = new StsComboBoxFieldBean(this, "edit", "Edit actions", editStrings);

    StsWell well = null;
    public String edit = NO_ACTION;

    static final String NO_ACTION = "No action";
    static final String NEW_LOG = "New log";
    static final String ADD_POINTS = "Add point";
    static final String APPLY_POINTS = "Apply";
    static final String CANCEL_POINTS = "Cancel";
    static final String[] editStrings = new String[] { NO_ACTION, NEW_LOG, ADD_POINTS, APPLY_POINTS, CANCEL_POINTS };

    static final StsWell nullWell = StsWell.nullWellConstructor("none");
    static final StsLogCurve nullLogCurve = StsLogCurve.nullLogCurveConstructor("none");

//    public static int displayableHeight = 600;

    public StsEditLogCurvePanel(StsEditLogCurve editLogCurve)
    {
        this.editLogCurve = editLogCurve;
        wizard = editLogCurve.getWizard();
        initializeWellPanel();
		//wellPanel.setPreferredSize(StsWellViewModel.getDefaultWellWindowSize());
    }

	public void initialize(Object[] wells)
	 {
		 wellListBean.setListItems(wells);
//		 logCurveListBean.setListItems(new Object[] {StsLogCurve.getNullLogCurve()});
	 }

	public void setSelectedWellInList(Object wellObject)
	{
        if(wellObject == null) return;
		if(wellListBean == null) return;
        wellListBean.setSelectedItem(wellObject);
	}

    public void LogCurveAdded(StsLogCurve logCurve)
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
        wellPanel.gbc.weighty = 0.0;
	    wellPanel.addEndRow(wellListBean);
        wellPanel.addEndRow(logCurveListBean);
        wellPanel.addEndRow(editBean);
        wellPanel.setBorder(BorderFactory.createEtchedBorder());
    }


    public void reshape(int x, int y, int width, int height)
    {
        super.reshape(x, y, width, height);
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
            if(wellObject == well) return;
            setSelectedWellInList(wellObject);
            // logCurve.setSelectedWell(wellObject);
            well = (StsWell) wellObject;
//System.out.println("well is "+well.getName());
            if (wellViewModel != null)
            {
                wellPanel.remove(wellViewModel.wellWindowPanel);
                wellViewModel = null;
            }
            if (actionManager != null)
			    actionManager.endCurrentAction();
			if (!well.getName().equals("none"))
			{
				wellViewModel = new StsWellWizardViewModel(well, wizard);
				System.out.println("well view model = " + wellViewModel + " " + well.wellViewModel);
				// setToPreviousWellViewModel();
				wellPanel.gbc.fill = GridBagConstraints.VERTICAL;
                wellPanel.gbc.weighty = 1.0;
				//wellPanel.gbc.gridwidth = 2;   // Combobox requires 2 grids width
				wellPanel.addEndRow(wellViewModel.wellWindowPanel);
				wellViewModel.initAndStartGL();

				StsObjectRefList logCurvesList = well.getLogCurves();
				logCurveListBean.setListItems(logCurvesList);
                checkChangeLogCurve(well);
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

    /** We have a new well.  If there was an old well, and we were editing a log of a particular type,
     *  switch to this same type in the new well.  Otherwise initialize to the first logCurve in the new well.
     * @param newWell  the newWell whose curve(s) we are now editing
     * @return
     */
    private StsLogCurve checkChangeLogCurve(StsWell newWell)
    {
        StsObjectRefList logCurves = newWell.getLogCurves();
        StsLogCurve newLogCurve = (StsLogCurve)logCurves.getObjectWithName(logCurve);
        if(newLogCurve != null)
            logCurve = newLogCurve;
        else
            logCurve = (StsLogCurve)logCurves.getFirst();
        if(!wellViewModel.isDisplayingLogCurve(logCurve))
        {
            wellViewModel.addLogCurve(logCurve);
        }
        return logCurve;
    }

    public Object getSelectedWell()
    {
        return well;
    }

    /** This log track might just have the single log being edited in it; we are now changing to a new log.
     *  If we delete the old log, the logTrack will be destroyed and must be rebuilt.
     *  Smoother transition is to add the new log and then delete the old log.
     * @param logCurve
     */
    public void setLogCurve(Object logCurve)
	{
        if(wellViewModel == null) return;
        if(this.logCurve == logCurve) return;
        StsLogCurve oldLogCurve = this.logCurve;
		this.logCurve = (StsLogCurve)logCurve;
        if(logCurve == nullLogCurve) return;
        if(this.logCurve != null)
        {
            this.logCurve.checkLoadVectors();
            StsLogCurvesView logCurvesView = wellViewModel.addCurveToLogTrack(this.logCurve);
            logCurvesView.getActionManager().addNewAction(editLogCurve);
        }
        if(oldLogCurve != null)
            wellViewModel.removeCurveFromLogTrack(oldLogCurve);
        rebuild();
        repaint();
    }

	public Object getLogCurve()
	{
	    return logCurve;
    }

    public void rebuild()
    {
        wizard.rebuild();
    }

    public String getEdit()
    {
        return edit;
    }

    public void setEdit(String edit)
    {
        if(editLogCurve == null) return;
        this.edit = edit;
        editLogCurve.setEdit(edit, logCurve, well);
    }
}
