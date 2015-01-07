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
import com.Sts.Actions.Wizards.StsWizardStep;
import com.Sts.Actions.Wizards.WizardHeaders.StsHeaderPanel;
import com.Sts.DBTypes.*;
import com.Sts.Interfaces.StsDialogFace;
import com.Sts.MVC.StsGLPanel;
import com.Sts.MVC.StsMessageFiles;
import com.Sts.MVC.View3d.StsGLPanel3d;
import com.Sts.MVC.View3d.StsView;
import com.Sts.MVC.ViewWell.StsLogCurvesView;
import com.Sts.Types.StsMouse;
import com.Sts.Types.StsMousePoint;
import com.Sts.Types.StsPoint;
import com.Sts.UI.Beans.StsGroupBox;
import com.Sts.UI.Beans.StsJPanel;
import com.Sts.UI.Beans.StsStringFieldBean;
import com.Sts.UI.StsOkCancelDialog;
import com.Sts.Utilities.StsException;
import com.Sts.Utilities.StsJOGLPick;
import com.Sts.Utilities.StsPickItem;
import com.Sts.Utilities.StsToolkit;

import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;

public class StsEditLogCurve extends StsWizardStep
{
    ArrayList<StsPoint> drawPoints = new ArrayList<StsPoint>();
    String edit = NO_ACTION;

    static final String NO_ACTION = StsEditLogCurvePanel.NO_ACTION;
    static final String NEW_LOG = StsEditLogCurvePanel.NEW_LOG;
    static final String ADD_POINTS = StsEditLogCurvePanel.ADD_POINTS;
    static final String APPLY_POINTS = StsEditLogCurvePanel.APPLY_POINTS;
    static final String CANCEL_POINTS = StsEditLogCurvePanel.CANCEL_POINTS;

    public StsEditLogCurve()
    {
        try
        {
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    StsEditLogCurvePanel panel;
    StsHeaderPanel header, info;

    StsWell selectedWell = null;
    StsLogCurve logCurve = null;

    static final boolean debug = false;

    public StsEditLogCurve(StsWizard wizard)
    {
        super(wizard);
        panel = new StsEditLogCurvePanel(this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
//        int height = header.displayableHeight + panel.displayableHeight;
//        panel.setPreferredSize(new Dimension(300, height));
        header.setTitle("Log Curve Editing");
        header.setSubtitle("Selecting Well");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#EditTD");
        header.setInfoText(wizardDialog, "(1) Select the well from list or in 3D window to edit.\n" +
                "(2) Select the log you wish to edit.\n" +
                " (3)Click and drag point or interval to adjust curve. Intermediate values will be set accordingly.\n");
    }

    public boolean start()
    {
        panel.initialize(getVisibleWellList());
        wizard.enableFinish();
        StsPoint.setCompareIndex(StsPoint.mIndex);
        return true;
    }

    public boolean end()
    {
        if(logCurve == null) return true;
        //logCurve.endTransaction();
        //logCurve.deleteSeismicCurtain();
        //model.removeDisplayableInstance(logCurve);
        return true;
    }

    public Object[] getVisibleWellList()
    {
        Object[] wells = model.getVisibleObjectList(StsWell.class);
        int nWells = wells.length;
        Object[] wellList = new Object[nWells + 1];
        StsWell nullWell = StsWell.nullWellConstructor("none");
        wellList[0] = nullWell;
        for(int n = 0; n < nWells; n++)
        {
            wellList[n + 1] = wells[n];
        }
        setSelectedWell(wellList[0]);
        return wellList;
    }

    public void setSelectedWell(Object object)
    {
        StsWell newSelectedWell = (StsWell) object;
        //if(logCurve != null) logCurve.deleteSeismicCurtain();
        if(newSelectedWell.getName().equals("none"))
            selectedWell = null;
        else
            selectedWell = newSelectedWell;
        panel.setSelectedWell(object);
        model.win3dDisplay();
    }

    public Object getSelectedWell()
    {
        return selectedWell;
    }

    public void setLogCurve(Object logCurve)
    {
        if(logCurve == null) return;
        this.logCurve = (StsLogCurve) logCurve;
    }

    public Object getLogCurve()
    {
        return logCurve;
    }

    /**
     * mouse action for 3d window
     */

    public boolean performMouseAction(StsMouse mouse, StsView view)
    {
        try
        {
            if(!(view instanceof StsLogCurvesView)) return false;
            if(panel.edit != ADD_POINTS) return false;
            StsLogCurvesView logCurvesView = (StsLogCurvesView)view;

            int leftButtonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);
            StsMousePoint mousePointGL = view.glPanel.getMousePointGL();
            if(leftButtonState == StsMouse.RELEASED)
            {
                float pointMDepth = (float)panel.wellViewModel.getMdepthFromGLMouseY(mousePointGL.y);
                float value = (float)getValueFromPanelX(mousePointGL.x, view);
                drawPoints.add(new StsPoint(value, pointMDepth));
                logCurvesView.setDrawPoints(drawPoints);
                panel.wellViewModel.wellWindowPanel.repaint();
                return true;
            }
            return false;
        }
        catch(Exception e)
        {
            StsException.outputException("Failed during pick.", e,
                    StsException.WARNING);
            return false;
        }
    }

    private double getValueFromPanelX(float x, StsView view)
    {
        return panel.logCurve.getValueFromPanelXFraction(x/view.glPanel.getWidth());
    }

    protected void setEdit(String edit, StsLogCurve logCurve, StsWell well)
    {
        this.edit = edit;
        if(edit == NEW_LOG)
        {
            createNewLog();
        }
        if(edit == APPLY_POINTS)
        {
            logCurve.applyPoints(drawPoints, well);
            drawPoints.clear();
        }
        else if(edit == CANCEL_POINTS)
            drawPoints.clear();
    }

    private void createNewLog()
    {
        new StsOkCancelDialog(null, new NewLog(), "Create New Log", true);
    }

    class NewLog implements StsDialogFace
    {
        StsJPanel panel = new StsJPanel();
        String panelName;
        String logNameString;

        void NewLog(String panelName)
        {
            this.panelName = panelName;
            panel.setName(panelName);
        }

        public void dialogSelectionType(int type)
        {
            boolean okSelected = (type == StsDialogFace.OK);
           // if(!okSelected) string = initialString;

        }
        public Component getPanel(boolean val) { return getPanel(); }

        public Component getPanel()
        {
            StsStringFieldBean logTypeStringBean = new StsStringFieldBean(this, "logNameString", "Log name");
            StsGroupBox groupBox = new StsGroupBox(panelName);
            groupBox.add(logTypeStringBean);
            return groupBox;
        }

        public StsDialogFace getEditableCopy()
        {
            return (StsDialogFace) StsToolkit.copyObjectNonTransientFields(this);
        }

        public void setLogNameString(String s) { logNameString = s; }
        public String getString() { return logNameString; }
    }
}
