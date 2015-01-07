package com.Sts.Actions.Wizards.SurfacesFromMarkers;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;

import javax.swing.*;
import java.awt.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsSelectMarkersPanel extends StsJPanel
{
    private StsSurfacesFromMarkersWizard wizard;
    private StsSelectMarkers wizardStep;

    private StsModel model = null;
    private StsMarker[] selectedMarkers = null;

    StsJPanel parameterPanel = StsJPanel.addInsets();

    JList markerList = new JList();
    DefaultListModel markerListModel = new DefaultListModel();
    protected JScrollPane jScrollPane1;

    StsDoubleFieldBean xOriginBean;
    StsDoubleFieldBean yOriginBean;
    StsFloatFieldBean xIncBean;
    StsFloatFieldBean yIncBean;
	StsIntFieldBean nRowsBean;
	StsIntFieldBean nColsBean;

    StsMarker[] markers;
    GridBagLayout gridBagLayout1 = new GridBagLayout();

    public StsSelectMarkersPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsSurfacesFromMarkersWizard)wizard;
        this.wizardStep = (StsSelectMarkers)wizardStep;
        try
        {
            constructBeans();
            jbInit();
            initialize();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    private void initialize()
    {
        model = wizard.getModel();
        markers = (StsMarker[])model.getCastObjectList(StsMarker.class);
        int nMarkers = markers.length;
        for(int n = 0; n < nMarkers; n++)
            markerListModel.addElement(markers[n].getName());
        markerList.setModel(markerListModel);
		parameterPanel.setPanelObject(wizard);
    }

    public void constructBeans()
    {
        xOriginBean = new StsDoubleFieldBean(wizard, "xOrigin", true, "X Origin:", true);
        yOriginBean = new StsDoubleFieldBean(wizard, "yOrigin", true, "Y Origin:", true);
        xIncBean = new StsFloatFieldBean(wizard, "xInc", true, "X Interval (m/f):", true);
        yIncBean = new StsFloatFieldBean(wizard, "yInc", true, "Y Interval (m/f):", true);
		nRowsBean = new StsIntFieldBean(wizard, "nRows", true, "# of Y Cells:", true);
		nColsBean = new StsIntFieldBean(wizard, "nCols", true, "# of X Cells:", true);
    }

    public StsMarker[] getSelectedMarkers()
    {
        int cnt = 0;
        if(markerList.isSelectionEmpty()) return null;
        int[] selectedIndices = markerList.getSelectedIndices();
        StsMarker[] selectedMarkers = new StsMarker[selectedIndices.length];
        for(int i=0; i< selectedIndices.length; i++)
            selectedMarkers[cnt++] = markers[selectedIndices[i]];
        return selectedMarkers;
    }

    void jbInit() throws Exception
    {
        this.gbc.fill = gbc.BOTH;
        this.gbc.anchor = gbc.WEST;
        this.gbc.weighty = 1.0;

        markerList.setBorder(BorderFactory.createEtchedBorder());
        jScrollPane1 = new JScrollPane(markerList);
        jScrollPane1.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        jScrollPane1.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);
        addEndRow(jScrollPane1);

        this.gbc.fill = gbc.HORIZONTAL;
        this.gbc.weighty = 0.0;

        parameterPanel.gbc.fill = gbc.HORIZONTAL;
        parameterPanel.gbc.anchor = gbc.WEST;

		parameterPanel.addToRow(nColsBean);
		parameterPanel.addEndRow(nRowsBean);
        parameterPanel.addToRow(xOriginBean);
		parameterPanel.addEndRow(yOriginBean);
        parameterPanel.addToRow(xIncBean);
        parameterPanel.addEndRow(yIncBean);

        addEndRow(parameterPanel);
    }
}
