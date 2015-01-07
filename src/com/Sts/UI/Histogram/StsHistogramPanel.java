//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI.Histogram;

import com.Sts.DBTypes.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.text.*;

public class StsHistogramPanel extends StsJPanel implements ActionListener
{
    public static final int VERTICAL = 0;
    public static final int HORIZONTAL = 1;
    int orientation = VERTICAL;
    boolean addLabels = false;
    // total range of display
    private float dataMin = 0;
    private float dataMax = 100;
    // range of data which is clipped (outside of clip is isVisible grey)
    private float clipDataMin = 0;
    private float clipDataMax = 100;
    // percent of data clipped off of min and max ends of spectrum
    private float minClipPercent = 100;
    private float maxClipPercent = 100;
    // beans with min and max percent clips isVisible (non-editable)
    StsJPanel minPanel = new StsJPanel();
    StsJPanel maxPanel = new StsJPanel();
    StsFloatFieldBean dataMinBean;
    StsFloatFieldBean dataMaxBean;
    StsFloatFieldBean clipDataMinBean;
    StsFloatFieldBean clipDataMaxBean;
    StsFloatFieldBean minClipPercentBean;
    StsFloatFieldBean maxClipPercentBean;

    // maximum amplitude of histogram
    float maxAmplitude;
    double histogramSum;
    ActionListener[] actionListeners = null;

    float verticalScale = -1.0f;
    StsJPanel histogramPanel = new StsJPanel();
    DecimalFormat colorFormat = new DecimalFormat("#.0##");
    Insets insets = null;
    Font defaultFont = new Font("Dialog", 0, 11);
    Font floatFont = new Font("Dialog", 0, 9);
    FontMetrics fm = getFontMetrics(floatFont);
    int dataX = 0, dataY = 0, dataXOld = -999, dataYOld = -999, dataWidth = 40, dataHeight = 20;
    Image offScreenImage = null;
    Graphics offScreenGraphics;
    Float dataValue = null;
    Float dataPercent = null;
    float[] opacityValues = null;
    int[] opacityX = new int[257];
    int[] opacityY = new int[257];
    StsColorscale colorscale = null;
    StsColorscalePanel colorscalePanel = null;
    private DecimalFormat numberFormat = new DecimalFormat("##0.0#");
    private JMenuItem kc = new JMenuItem("Set to Fully Opaque");

    int mouseBtn = 0;
    int previousIdx = -1;

    int maxWidth = 50;
    int nPoints = 255;
    int heightPerInterval = 1;
    int barLength = 255;
    int verticalBarWidth = 50;
    int horizontalBarHeight = 60;
    Rectangle barRect;
    int xOrigin, yOrigin;
    int minVal = 1;
    int maxVal = 255;
    int minIdx = 0;
    int maxIdx = 254;
    float topClip = 0.0f;
    float btmClip = 0.0f;

    int numSegments = 4;
    boolean segmentsOn = false;

    float[] histData = null;
    float[] histogramSamples = null;

    public StsHistogramPanel()
    {
        this(VERTICAL, false);
    }

    public StsHistogramPanel(int orient)
    {
        this(orient, false);
    }

    public StsHistogramPanel(int orient, boolean addLabels)
    {
        orientation = orient;
        this.addLabels = addLabels;
        initialize();
    }

    public StsHistogramPanel(int orient, float[] histData, float dataMin, float dataMax, boolean addLabels)
    {
        this(orient, addLabels);
        initializeData(histData, dataMin, dataMax);
    }

    public void initialize()
    {
        try
        {
            if(addLabels)
                buildPanelWithLabels();
            else
                buildPanel();
            initializeListeners();
        }
        catch (Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void updateData(StsSeismicBoundingBox volume)
    {
        updateData(volume.dataHist, volume.dataMin, volume.dataMax);
    }
    public void updateData(float[] histData, float dataMin, float dataMax)
    {
        if(histData == null) return;
        if(histData.length != 255)
        {
            System.err.println("Data Histogram must have 255 elements");
            return;
        }
        initializeData(histData, dataMin, dataMax);
        updateBeanValues();
    }

    public void updateSamples(float[] histogramSamples)
    {
        this.histogramSamples = histogramSamples;
        if(histData == null) return;
        if(histData.length != 255)
        {
            System.err.println("Data Histogram must have 255 elements");
            return;
        }
        initializeData(histData, dataMin, dataMax);
        updateBeanValues();
    }

    private void updateBeanValues()
    {
        if(addLabels)
        {
            dataMinBean.setValue(dataMin);
            dataMaxBean.setValue(dataMax);
        }
        setClipDataMinValue(dataMin);
        setClipDataMaxValue(dataMax);
    }

    private void updateBeanValues(StsSeismicBoundingBox volume)
    {
        if(addLabels)
        {
            dataMinBean.setValue(volume.dataMin);
            dataMaxBean.setValue(volume.dataMax);
        }
        setClipDataMinValue(volume.dataMin);
        setClipDataMaxValue(volume.dataMax);
    }

    private void initializeData(float[] histData, float dataMin, float dataMax)
    {
        this.histData = histData;
        this.dataMin = dataMin;
        this.dataMax = dataMax;
        maxAmplitude = 0;
        histogramSum = 0.0;
        for (int i = 0; i < histData.length; i++)
        {
            float data = histData[i];
            histogramSum += data;
            if (data > maxAmplitude)
                maxAmplitude = data;
        }
        setVerticalScale(maxAmplitude);
        resetOpacityValues();
    }

   private void buildPanel() throws Exception
    {
        histogramPanel.setFont(floatFont);
        gbc.insets = new Insets(0, 2, 0, 2);
        if(orientation == HORIZONTAL)
        {
            setHistogramPanelSize(barLength, horizontalBarHeight);
            gbc.anchor = gbc.WEST;
        }
        else
        {
            setHistogramPanelSize(verticalBarWidth, barLength);
        }
        add(histogramPanel);
    /*
        this.add(histogramPanel, new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            , GridBagConstraints.WEST, GridBagConstraints.NONE,
            new Insets(2, 2, 2, 2), 0, 0));
    */
    }

    private void buildPanelWithLabels() throws Exception
    {
        histogramPanel.setFont(floatFont);
        buildBeans();
 //       minPanel.setPreferredSize(150, 100);
        minPanel.add(dataMinBean);
        minPanel.add(clipDataMinBean);
        minPanel.add(minClipPercentBean);
//        maxPanel.setPreferredSize(150, 100);
        maxPanel.add(dataMaxBean);
        maxPanel.add(clipDataMaxBean);
        maxPanel.add(maxClipPercentBean);
        if(orientation == HORIZONTAL)
        {
            setHistogramPanelSize(barLength, horizontalBarHeight);
//            gbc.anchor = gbc.WEST;
            gbc.fill = gbc.HORIZONTAL;
            addToRow(minPanel, 1, 1);
            addToRow(histogramPanel, 1, 0);
            addEndRow(maxPanel, 1, 1);
        }
        else
        {
            setHistogramPanelSize(verticalBarWidth, barLength);
            gbc.fill = gbc.HORIZONTAL;
            add(minPanel, 1, 1);
            gbc.fill = gbc.NONE;
            add(histogramPanel, 1, 0);
            gbc.fill = gbc.HORIZONTAL;
            add(maxPanel, 1, 1);
        }
    }

    public int getMinMaxPanelHeight() { return minPanel.getBounds().height; }
    public void setClipEditable(boolean clip)
    {
        clipDataMinBean.setEditable(clip);
        clipDataMaxBean.setEditable(clip);
        if(!clip)
        {
            clipDataMin = dataMin;
            clipDataMax = dataMax;
            clipDataMinBean.setValue(dataMin);
            clipDataMaxBean.setValue(dataMax);
        }
    }

    private void setHistogramPanelSize(int width, int height)
    {
        Dimension size = new Dimension(width, height);
//        histogramPanel.setSize(size);
        histogramPanel.setMaximumSize(size);
        histogramPanel.setMinimumSize(size);
        histogramPanel.setPreferredSize(size);
    }

    private void buildBeans()
    {
        dataMinBean = new StsFloatFieldBean(this, "dataMin", false, "Min:");
        dataMaxBean = new StsFloatFieldBean(this, "dataMax", false, "Max:");
        clipDataMinBean = new StsFloatFieldBean(this, "clipDataMin", "Clip Min:");
        clipDataMaxBean = new StsFloatFieldBean(this, "clipDataMax", "ClipMax:");
        minClipPercentBean = new StsFloatFieldBean(this, "minClipPercent", false, "Min Clip %:");
        maxClipPercentBean = new StsFloatFieldBean(this, "maxClipPercent", false, "Max Clip %:");
    }

    private void initializeListeners()
    {
        addMouseListener();
        addMouseMotionListener();
        kc.addActionListener(this);
    }

    public void setDataRange(float dataMin, float dataMax)
    {
        this.setDataMin(dataMin);
        this.setDataMax(dataMax);
        if(addLabels)
        {
            dataMinBean.setValue(dataMin);
            dataMaxBean.setValue(dataMax);
            clipDataMinBean.setValue(dataMin);
            clipDataMaxBean.setValue(dataMax);
        }
        invalidate();
        repaint();
    }

    public void setDisplayIndexRange(int min, int max)
    {
        minVal = min;
        maxVal = max;
    }

    /** Programmatically set clipDataMin value; requires programmatically setting bean value */
    public void setClipDataMinValue(float clipDataMin)
    {
        if( histData == null) return;
        if(clipDataMinBean != null) clipDataMinBean.setValue(clipDataMin);
        setClipDataMin(clipDataMin);
    }

    public void setClipDataMin(float clipDataMin)
    {
        if( histData == null) return;
        this.clipDataMin = clipDataMin;
        minIdx = (int)(255*(clipDataMin - getDataMin())/(getDataMax() - getDataMin()));
        minClipPercent = 0;
        for(int n = 0; n < minIdx; n++)
            minClipPercent += histData[n];
        minClipPercent *= 100/histogramSum;
        if(minClipPercentBean != null) minClipPercentBean.setValue(minClipPercent);
        invalidate();
        repaint();
    }

    /** Programmatically set clipDataMin value; requires programmatically setting bean value */
    public void setClipDataMaxValue(float clipDataMax)
    {
        if( histData == null) return;
        if(clipDataMaxBean != null) clipDataMaxBean.setValue(clipDataMax);
        setClipDataMax(clipDataMax);
    }

    public void setClipDataMax(float clipDataMax)
    {
        if( histData == null) return;
        this.clipDataMax = clipDataMax;
        maxIdx = (int)(255*(clipDataMax - getDataMin())/(getDataMax() - getDataMin()));
        maxClipPercent = 0;
        for(int n = 254; n >= maxIdx; n--)
            maxClipPercent += histData[n];
        maxClipPercent *= 100/histogramSum;
        if(maxClipPercentBean != null) maxClipPercentBean.setValue(maxClipPercent);
        invalidate();
        repaint();
    }

    public void setColorscalePanel(StsColorscalePanel stsColorscalePanel)
    {
        this.colorscale = stsColorscalePanel.getColorscale();
        this.colorscalePanel = stsColorscalePanel;
        resetOpacityValues();
    }

    public void resetOpacityValues()
    {
        barRect = histogramPanel.getBounds();

        opacityX[0] = barRect.x + 2;
        opacityY[0] = nPoints;
        barRect = histogramPanel.getBounds();
        int bwidth = barRect.width;
        if(bwidth == 0)
            bwidth = verticalBarWidth;

		if(colorscale == null)
		{
			if (opacityValues == null)opacityValues = new float[nPoints + 1];
			for (int i = 1; i < nPoints + 1; i++)
				opacityValues[i - 1] = 1.0f;
		}
		else
        {
            if (opacityValues == null)opacityValues = new float[nPoints + 1];
            int[] colorscaleOpacityValues = colorscale.getOpacityValues();
            for (int i = 1; i < nPoints; i++)
				opacityValues[i] = colorscaleOpacityValues[i]/255.0f;
        }

        for (int i = 1; i < nPoints + 1; i++)
		{
			opacityX[i] = (int)(bwidth * opacityValues[i-1]) + 2;
            opacityY[i] = nPoints - i + 1;
        }
        opacityX[nPoints + 1] = barRect.x + 2;
        opacityY[nPoints + 1] = 0;
        repaint();
    }

    public void clearOpacityValues()
    {
        barRect = histogramPanel.getBounds();

        opacityX[0] = barRect.x + 1;
        opacityY[0] = nPoints + 1;
        barRect = histogramPanel.getBounds();
        int bwidth = barRect.width;
        if(bwidth == 0)
            bwidth = verticalBarWidth;

        for (int i = colorscalePanel.minSliderIndex + 1; i < colorscalePanel.maxSliderIndex + 1; i++)
        {
            if(colorscale == null)
                opacityValues[i - 1] = 1.0f;
            else
            {
                Color[] colors = colorscale.getNewColors();
                opacityValues[i - 1] = 1.0f;
            }
            opacityX[i] = (int)(bwidth * opacityValues[i-1]) + 1;
            opacityY[i] = nPoints - i + 1;
        }
        opacityX[nPoints + 1] = barRect.x + 1;
        opacityY[nPoints + 1] = 0;
        repaint();

        if (colorscale != null)
            colorscalePanel.setOpacityValues(opacityValues);
    }

    public float getMinClipPercent()
    {
        return minClipPercent;
    }

    public float getMaxClipPercent()
    {
        return maxClipPercent;
    }

    public void setVerticalScale(float scale)
    {
        verticalScale = scale;
    }

    public void actionPerformed(ActionEvent e)
    {
        int i;
        Object source = e.getSource();
        // Process Radio Button Input
        if(source == kc)
            if(colorscale != null) clearOpacityValues();
    }
/*
    public void repaint()
    {
        paint(getGraphics());
    }
*/
    public void paint(Graphics g)
    {
        if (g == null) return;
        if (histData == null) histData = new float[255];

        super.paint(g);

        barRect = histogramPanel.getBounds();
        if(barRect.width == 0 || barRect.height == 0) return;
        xOrigin = barRect.x;
        yOrigin = 0;

        Rectangle back = this.getVisibleRect();
        g.setColor(getBackground());
        g.fillRect(back.x, back.y, back.width, back.height);
        paintChildren(g);

        int bwidth = barRect.width;
        if (offScreenImage == null)
        {
            dataWidth = fm.stringWidth("99999.99");
            dataWidth += dataWidth * 0.05;
            if (dataWidth > (bwidth - 2))
            {
                dataWidth = bwidth - 2;
            }
//            dataHeight = fm.getHeight() * 2;
            dataHeight = fm.getHeight();
            dataHeight += dataHeight * 0.05;
            offScreenImage = createImage(dataWidth, dataHeight);
        }
        offScreenGraphics = offScreenImage.getGraphics();

        // Draw the histogram
        int x = xOrigin;
        int y = yOrigin;
        if (orientation == HORIZONTAL)
        {
            x = xOrigin;
            y = yOrigin;
            bwidth = barRect.height;
        }

        // Draw Color Scale - Bottom up
        float scale = bwidth / getClipDataMax();
        topClip = 0.0f;
        btmClip = 0.0f;
        int x1 = 0;
        int y1 = 0;
        int width = 0;
        int length = 0;
        int idx = 0;

        if(segmentsOn)
        {
            int segmentSize = heightPerInterval * histData.length / numSegments + 1;
            for(int i=0; i<numSegments; i++)
            {
                float colorNum = 1.0f - (float)i * 0.1f;
                g.setColor(new Color(colorNum,colorNum,colorNum));
                if(orientation == HORIZONTAL)
                {
                    g.fillRect(x + (i*segmentSize), y,  segmentSize, bwidth);
                }
                else
                {
                    g.fillRect(x, y + (i*segmentSize), bwidth, segmentSize);
                }
            }
        }

        if(verticalScale != -1.0f)
        {
            scale = bwidth / verticalScale;
        }

        for (int i = histData.length - 1; i >= 0; i--)
        {
            if (orientation == VERTICAL)
            {
                idx = 254 - i;
            }
            else
            {
                idx = i;

            }
            if(verticalScale != -1.0f)
            {
                if(histData[idx] > verticalScale)
                    histData[idx] = verticalScale;
            }
            if (idx > maxIdx)
            {
                g.setColor(Color.GRAY);
                topClip += histData[idx];
            }
            else if (idx < minIdx)
            {
                g.setColor(Color.GRAY);
                btmClip += histData[idx];
            }
            else
            {
                g.setColor(Color.BLUE);

                // draw a rectangle for each bar
            }
            if (orientation == VERTICAL)
            {
                x1 = x + bwidth - (int) (histData[idx] * scale);
                y1 = y + (i * heightPerInterval * (nPoints / histData.length));
                width = (int) (histData[idx] * scale);
                length = (heightPerInterval * (nPoints / histData.length));
            }
            else
            {
                x1 = x + (i * heightPerInterval * (nPoints / histData.length));
                y1 = y + bwidth - (int) (histData[idx] * scale);
                width = (heightPerInterval * (nPoints / histData.length));
                length = (int) (histData[idx] * scale);
            }
            g.fillRect(x1, y1, width, length);
        }
        g.setColor(Color.black);
        if (orientation == HORIZONTAL)
        {
            g.drawRect(x, y, heightPerInterval * histData.length, bwidth);
        }
        else
        {
            g.drawRect(x, y, bwidth, heightPerInterval * histData.length);

        }
        if (mouseBtn == 1)
        {
            paintOffScreen(offScreenGraphics);
            g.drawImage(offScreenImage, dataX, dataY /*- fm.getHeight()*/, this);
        }
        if (orientation != HORIZONTAL && colorscale != null)
        {
            g.setColor(new Color(1.0f, 0.0f, 0.0f, 0.2f));
            g.fillPolygon(opacityX, opacityY, nPoints + 2);
            g.setColor(new Color(0.5f, 0.5f, 0.5f));
            g.drawPolyline(opacityX, opacityY, nPoints + 2);
        }
    }

    private void addKeyPopup(int x, int y)
    {
        JPopupMenu tp = new JPopupMenu("Histogram Popup");
        this.add(tp);
        tp.add(kc);
        tp.show(this, x, y);
    }

    private void addMouseListener()
    {
        MouseListener mouseListener = new MouseAdapter()
        {
            public void mouseReleased(MouseEvent e)
            {
                dataValue = null;
                dataPercent = null;
                if (mouseBtn == 3)
                {
                    if (colorscale != null)
                    {
                        colorscalePanel.setOpacityValues(opacityValues);
                    }
                    mouseBtn = 0;
                    previousIdx = -1;
                    repaint();
                }
            }
            public void mousePressed(MouseEvent e)
            {
                mouseBtn = e.getButton();
                int mods = e.getModifiers();
                if((e.isShiftDown() && ((mods & InputEvent.BUTTON3_MASK)) != 0) || ((mods & InputEvent.BUTTON2_MASK) != 0))
                {
                    if(colorscale != null)
                        addKeyPopup(e.getX(),e.getY());
                }
            }
        };
        this.addMouseListener(mouseListener);
    }

    private void addMouseMotionListener()
    {
        MouseMotionListener motionListener = new MouseMotionAdapter()
        {
            public void mouseDragged(MouseEvent e)
            {
                int y = 0;
                int pickedY = e.getY();
                int pickedX = e.getX();

                if (orientation == HORIZONTAL)
                {
                    if (pickedY > barRect.y + barRect.height)
                        pickedY = barRect.y + barRect.height;
                    y = - (pickedX - xOrigin);
                }
                else
                {
                    if (pickedX > barRect.x + barRect.width)
                        pickedX = barRect.x + barRect.width;
                    y = pickedY - (barRect.y + barRect.height);

                }
                int index = -y / heightPerInterval;

                if (mouseBtn == e.BUTTON1)
                {
                    // Output the data value associated with the mouse position. SAJ
                    if ( (index < 0) || (index > maxIdx))
                    {
                        dataValue = null;
                        dataPercent = null;
                        return;
                    }
                    else
                    {
                        offScreenGraphics.dispose();
                        dataPercent = new Float(colorFormat.format(histData[index]));
                        dataValue = new Float(colorFormat.format(getDataValue(index)));

                        dataX = e.getX();
                        if ( (dataWidth + dataX) > barRect.width)
                        {
                            dataX = barRect.width - dataWidth;
                        }
                        else if (dataX < barRect.x)
                        {
                            dataX = barRect.x + 1;

                        }
                        dataY = e.getY() - fm.getHeight();
                        if ( (e.getY() - fm.getHeight()) < barRect.y)
                        {
                            dataY = barRect.y + 3;
                        }
                        else if (e.getY() > barRect.y + barRect.height)
                        {
                            dataY = barRect.y + barRect.height - fm.getHeight() - 3;
                        }
                    }
                    repaint();
                }
                else if ( (mouseBtn == e.BUTTON3) && (orientation != HORIZONTAL))
                {
                    if((index < minIdx) || (index > maxIdx))
                        return;

                    dataX = e.getX() + barRect.x;
                    if (dataX > barRect.width + barRect.x)
                        dataX = barRect.width + barRect.x;
                    else if (dataX < barRect.x)
                        dataX = barRect.x;

                    dataY = e.getY();
                    if (dataY < barRect.y)
                        dataY = barRect.y;
                    else if (dataY > barRect.y + barRect.height)
                        dataY = barRect.y + barRect.height;

                    opacityX[index + 1] = dataX;
                    opacityY[index + 1] = dataY;
                    opacityValues[index] = (float) (dataX - barRect.x) / verticalBarWidth;
                    if(opacityValues[index] > 1.0) opacityValues[index] = 1.0f;
                    if(opacityValues[index] < 0.0f) opacityValues[index] = 0.0f;

                        // Interpolated between previous index and current one
                    if((previousIdx != -1) && (Math.abs(previousIdx - index) > 1))
                    {
                        int length = previousIdx - index;
                        float scale = (float)((float)(opacityX[index+1] - opacityX[previousIdx + 1])/Math.abs(length));
                        if(length > 0)
                            for(int n = 1; n < length; n++)
                            {
                                opacityX[index + n + 1] = opacityX[index + 1] + (int)(scale * (float)n);
                                opacityY[index + n + 1] = opacityY[index + 1] + -(heightPerInterval * n);
                                opacityValues[index + n] = (float)(((float)opacityX[index + 1] + (scale * (float)n) - barRect.x) / verticalBarWidth);
                                if(opacityValues[index+n] > 1.0) opacityValues[index+n] = 1.0f;
                                if(opacityValues[index+n] < 0.0f) opacityValues[index+n] = 0.0f;
                            }
                        else
                            for(int n = -1; n > length; n--)
                            {
                                opacityX[index + n + 1] = opacityX[index + 1] + (int)(scale * (float)n);
                                opacityY[index + n + 1] = opacityY[index + 1] + -(heightPerInterval * n);
                                opacityValues[index + n] = (float)(((float)opacityX[index + 1] + (scale * (float)n) - barRect.x) / verticalBarWidth);
                                if(opacityValues[index+n] > 1.0) opacityValues[index+n] = 1.0f;
                                if(opacityValues[index+n] < 0.0f) opacityValues[index+n] = 0.0f;
                            }
                    }
                    previousIdx = index;
                    repaint();
                }
            }
        };
        this.addMouseMotionListener(motionListener);
    }

        public float getDataValue(int idx)
        {
            float scale = (maxVal - minVal) / nPoints;
            float value = idx * scale + minVal;
            return value;
        }

        public int getNumberOfSegments() { return numSegments; }
        public void setNumberOfSegment(int num)
        {
            numSegments = num;
        }
        public boolean getSegmentsOn() { return segmentsOn; }
        public void setSegmentsOn(boolean bool)
        {
            segmentsOn = bool;
        }

        private void paintOffScreen(Graphics og)
        {
            if(dataValue == null) return;
            og.clearRect(0, 0, dataWidth, dataHeight);
            og.setColor(Color.red);
            og.drawRect(0, 0, dataWidth - 1, dataHeight - 1);
            og.setColor(Color.black);
            og.setFont(floatFont);
            int x = (int) ( (float) dataWidth * 0.05f);
//            int y = (int) ((float) dataHeight * 0.10f);
            int y = (int) ( (float) dataHeight - (float) dataHeight * 0.10f);
            og.drawString(dataPercent.toString(), x, y);
//            og.drawString(dataValue.toString(), x, y);
//            og.drawString(dataPercent.toString(), x, y - fm.getHeight());
        }

        public boolean setClip(int min, int max)
        {
            boolean changed = false;
            if (minIdx != min)
            {
                minIdx = min;
                changed = true;
            }
            if (maxIdx != max)
            {
                maxIdx = max;
                changed = true;
            }
            return changed;
        }

        public void clearAll()
        {
            minIdx = 0;
            maxIdx = 254;
            setClipDataMax(-1000.0f);
            setClipDataMin(1000.0f);
            topClip = 0.0f;
            btmClip = 0.0f;

            histData = null;
        }

        public int getNumberIndices()
        {
            return barLength;
        }

        public float[] getOpacityValues()
        {
            return opacityValues;
        }

        public float getTopPercentageClipped()
        {
            return topClip;
        }
        public float getBottomPercentageClipped()
        {
            return btmClip;
        }
        public float[] getDataHistogram()
        {
            return histData;
        }
        public int getHeightPerInterval()
        {return heightPerInterval;
        }
        public void setHeightPerInterval(int heightPerInterval)
        {this.heightPerInterval = heightPerInterval;
        }
        public int getVerticalBarWidth()
        {
            return verticalBarWidth;
        }
        public void setVerticalBarWidth(int verticalBarWidth)
        {
            this.verticalBarWidth = verticalBarWidth;
        }
        public int getNPoints()
        {return nPoints;
        }
        public void setNPoints(int nPoints)
        {this.nPoints = nPoints;
        }

        static public void main(String[] args)
        {
            float[] dataHist =
                {
                0.025f, 0.035f, 0.045f, 0.09f, 0.11f, 0.2f, 0.35f, 0.56f, 0.78f, 1.0f, 1.25f, 1.48f,
                1.67f, 2.34f, 3.56f, 3.24f, 3.45f, 3.12f, 2.76f, 2.34f, 2.1f, 1.5f, 1.25f, 1.03f,
                0.99f, 0.92f, 0.81f, 0.71f, 0.61f, 0.5f, 0.35f, 0.56f, 0.78f, 1.0f, 1.25f, 1.48f,
                1.67f, 2.34f, 3.56f, 3.24f, 3.45f, 3.12f, 2.76f, 2.34f, 2.1f, 1.5f, 1.25f, 1.03f,
                1.67f, 2.34f, 3.56f, 3.24f, 3.45f, 3.12f, 2.76f, 2.34f, 2.1f, 1.5f, 1.25f, 1.03f,
                0.025f, 0.035f, 0.045f, 0.09f, 0.11f, 0.2f, 0.35f, 0.56f, 0.78f, 1.0f, 1.25f, 1.48f,
                0.025f, 0.035f, 0.045f, 0.09f, 0.11f, 0.2f, 0.35f, 0.56f, 0.78f, 1.0f, 1.25f, 1.48f,
                1.67f, 2.34f, 3.56f, 3.24f, 3.45f, 3.12f, 2.76f, 2.34f, 2.1f, 1.5f, 1.25f, 1.03f,
                0.025f, 0.035f, 0.045f, 0.09f, 0.11f, 0.2f, 0.35f, 0.56f, 0.78f, 1.0f, 1.25f, 1.48f,
                1.67f, 2.34f, 3.56f, 3.24f, 3.45f, 3.12f, 2.76f, 2.34f, 2.1f, 1.5f, 1.25f, 1.03f,
                1.67f, 2.34f, 3.56f, 3.24f, 3.45f, 3.12f, 2.76f, 2.34f, 2.1f, 1.5f, 1.25f, 1.03f,
                0.025f, 0.035f, 0.045f, 0.09f, 0.11f, 0.2f, 0.35f, 0.56f, 0.78f, 1.0f, 1.25f, 1.48f,
                0.025f, 0.035f, 0.045f, 0.09f, 0.11f, 0.2f, 0.35f, 0.56f, 0.78f, 1.0f, 1.25f, 1.48f,
                1.67f, 2.34f, 3.56f, 3.24f, 3.45f, 3.12f, 2.76f, 2.34f, 2.1f, 1.5f, 1.25f, 1.03f,
                0.025f, 0.035f, 0.045f, 0.09f, 0.11f, 0.2f, 0.35f, 0.56f, 0.78f, 1.0f, 1.25f, 1.48f,
                1.67f, 2.34f, 3.56f, 3.24f, 3.45f, 3.12f, 2.76f, 2.34f, 2.1f, 1.5f, 1.25f, 1.03f,
                1.67f, 2.34f, 3.56f, 3.24f, 3.45f, 3.12f, 2.76f, 2.34f, 2.1f, 1.5f, 1.25f, 1.03f,
                0.025f, 0.035f, 0.045f, 0.09f, 0.11f, 0.2f, 0.35f, 0.56f, 0.78f, 1.0f, 1.25f, 1.48f,
                1.67f, 2.34f, 3.56f, 3.24f, 3.45f, 3.12f, 2.76f, 2.34f, 2.1f, 1.5f, 1.25f, 1.03f,
                0.025f, 0.035f, 0.045f, 0.09f, 0.11f, 0.2f, 0.35f, 0.56f, 0.78f, 1.0f, 1.25f, 1.48f,
                .9f, 0.75f, 0.23f, 0.04f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
                0.0f
            };

            try
            {
                StsHistogramPanel histogramPanel = new StsHistogramPanel(HORIZONTAL, dataHist, 0, 100, true);
                StsToolkit.createDialog(histogramPanel, false);

                StsToolkit.sleep(500);
                histogramPanel.setClipDataMinValue(25);
                StsToolkit.sleep(500);
                histogramPanel.setClipDataMaxValue(75);
            }
            catch (Exception e)
            {
                e.printStackTrace();
            }
        }

    public float getDataMin()
    {
        return dataMin;
    }

    public void setDataMin(float dataMin)
    {
        this.dataMin = dataMin;
    }

    public float getDataMax()
    {
        return dataMax;
    }

    public void setDataMax(float dataMax)
    {
        this.dataMax = dataMax;
    }

    public float getClipDataMin()
    {
        return clipDataMin;
    }

    public float getClipDataMax()
    {
        return clipDataMax;
    }
}
