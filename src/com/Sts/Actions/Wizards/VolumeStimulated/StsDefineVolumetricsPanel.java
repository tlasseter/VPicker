package com.Sts.Actions.Wizards.VolumeStimulated;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import info.monitorenter.gui.chart.*;
import info.monitorenter.gui.chart.traces.*;
import info.monitorenter.gui.chart.views.*;

import java.awt.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsDefineVolumetricsPanel extends StsJPanel
{
    private StsVolumeStimulatedWizard wizard;
    private StsDefineVolumetrics wizardStep;
 
    private StsModel model = null;

    StsGroupBox criteriaBox = new StsGroupBox("Define Rock Properties");
    StsFloatFieldBean porosityBean = new StsFloatFieldBean();
    StsFloatFieldBean permBean = new StsFloatFieldBean();
    StsFloatFieldBean areaScaleBean = new StsFloatFieldBean();
    
    StsGroupBox plotBox = new StsGroupBox("Estimated Production");
    StsJPanel plotPanel = new StsJPanel();
    Chart2D chart = null;
    ChartPanel chartPanel = null;
    ITrace2D trace = null;
    
	StsFractureSet[] selectedFractureSets = null;
	
    public StsDefineVolumetricsPanel(StsVolumeStimulatedWizard wizard, StsDefineVolumetrics wizardStep)
    {
        this.wizard = wizard;
        this.wizardStep = wizardStep;
        buildPanel();
    }

    public void buildPanel()
    {
        porosityBean.initialize(wizard, "porosity", 0.1f, 20.0f, "Effective Porosity:", true);
        porosityBean.fixStep(0.1f);
        porosityBean.setToolTipText("Specify the effective porosity of the combined natural fractures and rock matrix.");
        
        permBean.initialize(wizard, "permeability", 0.0f, 1.0f, "Effective Permeability:", true);
        permBean.fixStep(0.001f);
        permBean.setToolTipText("Specify the effective permeability of the combined natural fractures and rock matrix.");
        
        areaScaleBean.initialize(wizard, "areaScale", 1.0f, 25.0f, "Area Scalar:", true);
        areaScaleBean.fixStep(0.1f);
        areaScaleBean.setToolTipText("Specify the scalar to apply to the fracture area.");
        
        criteriaBox.gbc.fill = gbc.HORIZONTAL;
        criteriaBox.add(porosityBean);
        criteriaBox.add(permBean);
        criteriaBox.add(areaScaleBean);
        gbc.fill = gbc.HORIZONTAL;
        gbc.anchor = gbc.NORTH;
        gbc.weighty = 0.1;
        add(criteriaBox);
        
        plotPanel.setBackground(Color.WHITE);
        plotBox.gbc.fill = gbc.BOTH;
        plotBox.add(plotPanel);
        gbc.fill = gbc.BOTH;
        gbc.anchor = gbc.SOUTH;
        gbc.weighty = 1.0;       
        add(plotBox);
    }

    public void initialize()
    {
        model = wizard.getModel();
        gbc.fill = gbc.HORIZONTAL;

        selectedFractureSets = ((StsVolumeStimulatedWizard)wizard).getFractureSets();
        initializePlot();
    }
    
    public void initializePlot()
    {
        trace = new Trace2DSimple();

        chart = createChart();
        chart.setMinPaintLatency(100);

        //new LayoutFactory.BasicPropertyAdaptSupport(this.getRootPane(), chart);
        //LayoutFactory lf = LayoutFactory.getInstance();
        //lf.setShowTraceZindexMenu(false);
        //lf.setShowTraceNameMenu(false);
        //lf.setShowAxisXRangePolicyMenu(false);
        chart.setEnablePopup(false);
        chartPanel = new ChartPanel(chart);
        plotPanel.gbc.fill = gbc.BOTH;
        plotPanel.add(chartPanel);
    	chart.setRequestedRepaint(true);
        chart.updateUI();
        wizard.rebuild();
    }
    
    private Chart2D createChart()
    {
        Chart2D chart1 = new Chart2D();
        buildTraces(chart1);
        chart1.setPaintLabels(true);
        chart1.setPaintCrosshair(true);
        return chart1;
    }

    private void buildTraces(Chart2D chart1)
    {
        Color color = StsColor.BLACK.getColor();
        trace = createDataset();
        trace.setColor(color);
        trace.setVisible(true);
        chart1.addTrace(trace);
    }
    
    private Trace2DSimple createDataset()
    {
        Trace2DSimple trace = new Trace2DSimple("Estimate");
        for (int i = 0; i < 10; i++)
        {
            trace.addPoint(i, Math.cos(i));
        }
        return trace;
    }
}
