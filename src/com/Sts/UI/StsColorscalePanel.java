
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Histogram.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;

/** This is the UI for an StsColorscale and has the button and slider controls.
 *
 *  A view-only StsColorscalePanel can be added to a StsColorscaleFieldBean (not editable) for inclusion on the object panel.
 *  Or an editable StsColorscalePanel can be added to an StsEditableColorscalePanel which
 *  includes an edit button for inclusion on the object Panel.
 *
 *  An editable StsColorscalePanel can be inserted into an StsColorscaleDialog for an editable dialog UI.
 */

 public class StsColorscalePanel extends JPanel implements ActionListener, FocusListener, ChangeListener
{
    ActionListener[] actionListeners = null;
    int orientation = StsParameters.VERTICAL;  // if horizontal, buttons are to the left of vertical bar
    JPanel colorbarPanel = new JPanel();
    StsHistogramPanel histogramPanel = new StsHistogramPanel();
    JRadioButton transparencyBtn = new JRadioButton();
    JRadioButton opaqueBtn = new JRadioButton();
    JRadioButton clipBtn = new JRadioButton();
    JRadioButton compressBtn = new JRadioButton();
    JCheckBox flipBtn = new JCheckBox("Flip");
    JCheckBox syncBtn = new JCheckBox("Sync Sliders");
    JSlider rotateSlider = new JSlider(JSlider.HORIZONTAL);

    JCheckBox voxelBtn = new JCheckBox("Voxel");
	StsFloatFieldBean voxelMinBean = new StsFloatFieldBean(this, "voxelMin", true, null, false);
	StsFloatFieldBean voxelMaxBean = new StsFloatFieldBean(this, "voxelMax", true, null, false);
//    JTextField voxelMinText = new JTextField();
//    JTextField voxelMaxText = new JTextField();
    JButton voxelSetBtn = new JButton("Set");
    StsColorscale.VoxelKey currentVoxelKey = null;
    boolean showVoxelControls = true;

    float dataValue = nullValue;
    int dataX=0, dataY=0, dataWidth = 40, dataHeight = 20;
    int voxelStart=0, voxelEnd=0;
    Image offScreenImage = null;
    Graphics offScreenGraphics;

	StsFloatFieldBean maxSliderBean = new StsFloatFieldBean(); //this, "maxSliderValue", true, "Max");
	StsFloatFieldBean minSliderBean = new StsFloatFieldBean(); //this, "minSliderValue", true, "Min");
//    JTextField maxData = new JTextField();
//    JTextField minData = new JTextField();
//    JLabel maxLabel = new JLabel("Data Maximum:");
//    JLabel minLabel = new JLabel("Data Minimum:");
    private ButtonGroup transparencyGroup = new ButtonGroup();
    private ButtonGroup modeGroup = new ButtonGroup();
//    boolean refresh = false;
    boolean viewOnly = false;
    boolean labelsOn = true;
    boolean syncSliders = false;
    float[] originalOpacity = null;

    StsSpectrum spectrum = null;
    StsFieldBean fieldBean = null;

	/** another colorscalePanel (on objectPanel which needs to update in sync with this one */
	StsColorscalePanel observerColorscalePanel = null;

    Polygon keyPolygon;
    Polygon voxelPolygon;
    Color key = Color.red;
    Color originalColor = Color.black;
    static int keyHeight = 6;
    static int keyWidth = 22;
    int keyIndex = -1;
    int numKeyIndices = 1;
    private JMenuItem colorChangeKey = new JMenuItem("Change Key Color");
    private JMenuItem colorResetKey = new JMenuItem("Reset Key Color");
    private JMenuItem increaseWidthKey = new JMenuItem("Increase Key Width");
    private JMenuItem decreaseWidthKey = new JMenuItem("Decrease Key Width");

	static final float nullValue = StsParameters.nullValue;

    private JMenuItem vd = new JMenuItem("Delete Voxel Key");

    static JColorChooser jColorChoose = null;

	static final boolean debug = false;

    public int minSliderIndex;
    public int maxSliderIndex;

    Insets insets = null;
    int nPoints = 255;
    int maxIndex = 254;
    int heightPerInterval = 1;
    int barHeight;
    int barWidth = 50;
    Rectangle cbarRect, hbarRect, back;
    Insets cInsets;
    int xOrigin, yOrigin;
    Polygon minSliderPolygon, maxSliderPolygon;
    int pickedX, pickedY;
    int pickedIndex = 0;
    boolean leftSide = false;
    int flipX = 1;
    float maxValue = 1.0f;
    float minValue = 0.0f;
    float maxSliderValue = 0.25f;
    float minSliderValue = 0.75f;
//    int voxelMaxValue = 255;
//    int voxelMinValue = 0;
    float scale;

    Color[] newColors = null;

    int mode = COLORSCALE;
    public final static int COLORSCALE = 0;
    public final static int SPECTRUM = 1;
    public final static int COLORS = 2;
    public final static int BEAN = 3;

    byte sliderSelected = NONE;
    static private byte NONE = 0;
    static private byte MIN_SLIDER = -1;
    static private byte MAX_SLIDER = 1;
    static private byte KEY_SLIDER = 2;
    static private byte VOXEL_SLIDER = 3;

    static int sliderHeight = 5;
    static int sliderWidth = 32;

    Color[] inColors = new Color[255];

    Point[] p;

    /*  This defines minSlider geometry; note y increases down in Java graphics
     *
     *   p1            p0
     *   ---------------=============   ^           ^          y   x ------->
     *   |             | intervalHeight |           |          |
     *   |             |                v           |          |
     *   |            / p4                     totalHeight     |
     *   |           /                              |          |
     *   |----------/                               v          V
     *   p2        p3
     *
     *   |<-- width -->|
     *
     *   For maxslider, flip minSlider and 0,0 should be at p4, so include an intervalHeight offset.
     *   To put sliders on right-side, flip x values;
     */


    StsDecimalFormat labelFormat = new StsDecimalFormat(4);
    StsDecimalFormat colorFormat = new StsDecimalFormat(3);
    StsColorscale colorscale = null;

    int perColorHeight = 1;
    int maxWidth = 150;
    boolean hasHistogram = false;

    Font defaultFont = new Font("Dialog",0,11);
    Font floatFont = new Font("Dialog",0,9);
    FontMetrics fm = getFontMetrics(floatFont);

    // If this is a panel inside a fieldBean, get spectrum from fieldBean
    public StsColorscalePanel(StsFieldBean fieldBean)
    {
        this(fieldBean, StsParameters.VERTICAL);
    }

    public StsColorscalePanel(StsFieldBean fieldBean, int orientation)
    {
        this.viewOnly = true;
        this.fieldBean = fieldBean;
        this.orientation = orientation;

        this.spectrum = (StsSpectrum)fieldBean.getPanelObject();
        mode = BEAN;

        try
        {
            init();
            jbInit();
        }
        catch (Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public StsColorscalePanel()
    {
        viewOnly = false;

        try
        {
            init();
            jbInit();
        }
        catch (Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public StsColorscalePanel(boolean vo, int mode)
    {
        viewOnly = vo;
        hasHistogram = false;
        this.mode = mode;

        try
        {
            init();
            jbInit();
            if(!viewOnly) initButtons();
			if(colorscale != null) colorscale.setEditRange(minSliderValue, maxSliderValue);
            adjustColorscale();
        }
        catch (Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public StsColorscalePanel(boolean vo)
    {
        viewOnly = vo;

        try
        {
            init();
            jbInit();
            initButtons();
			if(colorscale != null) colorscale.setEditRange(minSliderValue, maxSliderValue);
            adjustColorscale();
        }
        catch (Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public StsColorscalePanel(Color[] colors)
    {
        this(true);
        this.inColors = colors;
        this.spectrum = null;
        this.mode = COLORS;
    }

    public StsColorscalePanel(StsSpectrum spectrum)
    {
        this(true);
        this.spectrum = spectrum;
        this.mode = SPECTRUM;
    }

    public StsColorscalePanel(StsColorscale colorscale)
    {
        this(colorscale, StsParameters.VERTICAL, false);
    }

	public StsColorscalePanel(StsColorscalePanel observerColorscalePanel)
	{
		this(observerColorscalePanel.getColorscale(), StsParameters.VERTICAL, false);
		this.observerColorscalePanel = observerColorscalePanel;
        float[] histData = observerColorscalePanel.getHistogram();
        float dataMin = observerColorscalePanel.minValue;
        float dataMax = observerColorscalePanel.maxValue;
        this.histogramPanel.updateData(histData, dataMin, dataMax);
    }

	public StsColorscalePanel(StsColorscale colorscale, boolean viewOnly)
	{
		this(colorscale, StsParameters.VERTICAL, viewOnly);
    }

	public StsColorscalePanel(StsColorscale colorscale, int orientation)
	{
		this(colorscale, orientation, false);
	}

    public StsColorscalePanel(StsColorscale colorscale, int orientation, boolean viewOnly)
    {
        this(colorscale, orientation, viewOnly, true);
    }

	public StsColorscalePanel(StsColorscale colorscale, int orientation, boolean viewOnly, boolean voxelControls)
    {
        this.colorscale = colorscale;
        this.orientation = orientation;
        this.showVoxelControls = voxelControls;

        this.spectrum = colorscale.getSpectrum();
        this.mode = COLORSCALE;
        this.viewOnly = viewOnly;

        try
        {
            init();
            jbInit();
            initButtons();
        }
        catch (Exception ex)
        {
            ex.printStackTrace();
        }
    }

    private void init()
    {
        if(!viewOnly)
        {
            addMouseListener();
            addMouseMotionListener();
        }

		initializeSliderScale();

        if(histogramPanel != null)
        {
            histogramPanel.setDataRange(minValue, maxValue);
            if(colorscale != null)
			{
				histogramPanel.setColorscalePanel(this);
// SAJ - This causes the opacity to be reset to 1.0 when an object is selected on the object panel
// not sure why it is here, commenting out fixed reset problem and don't see any adverse effects
//                histogramPanel.clearOpacityValues();
			}
        }
//        minData.setText(labelFormat.format(minValue));
//        maxData.setText(labelFormat.format(maxValue));

        initializeSliderPolygons();
        barHeight = (maxIndex+1)*heightPerInterval;

        if(colorscale != null)
        {
            newColors = colorscale.getNewColors();
            if(newColors == null) return;

            keyIndex = newColors.length / 2;
            originalColor = newColors[keyIndex];
            key = newColors[keyIndex];
            numKeyIndices = 0;

            float[] sliderRange = colorscale.getEditRange();
			checkSetSliderValues(sliderRange[0], sliderRange[1]);
//			float[] range = colorscale.range();
 //           minData.setText(Float.toString(range[0]));
 //           maxData.setText(Float.toString(range[1]));
 //           setMinSliderFromValue(range[0]);
 //           setMaxSliderFromValue(range[1]);

	        StsColorscale.VoxelKey[] voxelKeys = colorscale.getVoxelKeys();
			if (voxelKeys != null && voxelKeys.length > 0 && voxelKeys[0] != null)
			{
				currentVoxelKey = voxelKeys[0];
				setVoxelRange();
			}

            colorChangeKey.addActionListener(this);
            colorResetKey.addActionListener(this);
            increaseWidthKey.addActionListener(this);
            decreaseWidthKey.addActionListener(this);

            vd.addActionListener(this);

//            refresh = true;
			if(colorscale != null) colorscale.setEditRange(minSliderValue, maxSliderValue);
            initializeColorscale();
//            repaintModel();
        }
    }

    private void setVoxelRange()
    {
        if(currentVoxelKey == null)
        {
            float[] range = colorscale.getRange();
            float avg = (range[0] + range[1])/2;
			voxelMinBean.setValue(avg);
			voxelMaxBean.setValue(avg);
        }
        else
        {
            float[] voxelRange = currentVoxelKey.getDataRange();
			voxelMinBean.setValue(voxelRange[0]);
			voxelMaxBean.setValue(voxelRange[1]);
         }
    }

	public float getVoxelMin()
	{
		if(currentVoxelKey != null)
		{
			return currentVoxelKey.getDataRange()[0];
		}
		else
			return this.minSliderValue;
	}

	public float getVoxelMax()
	{
		if(currentVoxelKey != null)
		{
			return currentVoxelKey.getDataRange()[1];
		}
		else
			return this.maxSliderValue;
	}

	public void setVoxelMin(float value)
	{
		if(currentVoxelKey == null) return;
		currentVoxelKey.setMinValue(value);
		setVoxelTransparency();
	}

	public void setVoxelMax(float value)
	{
		if(currentVoxelKey == null) return;
		currentVoxelKey.setMaxValue(value);
		setVoxelTransparency();
	}

    private void initButtons()
    {
        if(colorscale != null)
        {
            flipBtn.setSelected(colorscale.getFlip());
            if(colorscale.getTransparencyMode())
                transparencyBtn.setSelected(true);
            else
                opaqueBtn.setSelected(true);
            if(colorscale.getCompressionMode() == colorscale.CLIPPED)
                clipBtn.setSelected(true);
            else
                compressBtn.setSelected(true);

			syncSliders = colorscale.getIsSynched();
			syncBtn.setSelected(syncSliders);

            voxelBtn.setSelected(false);
            voxelMinBean.setEnabled(false);
            voxelMaxBean.setEnabled(false);
            voxelSetBtn.setEnabled(false);
        }
        else
        {
            flipBtn.setSelected(false);
			syncBtn.setSelected(false);
            voxelBtn.setSelected(false);
            compressBtn.setSelected(true);
            transparencyBtn.setSelected(true);
        }
    }

    public void reinitialize()
    {
		float[] range = colorscale.getRange();
		minValue = range[0];
		maxValue = range[1];
		scale = (maxValue - minValue)/maxIndex;
		minSliderBean.setRange(minValue, maxValue);
		maxSliderBean.setRange(minValue, maxValue);
	    float[] editRange = colorscale.getEditRange();
		checkSetSliderValues(editRange[0], editRange[1]);
         initButtons();
		 adjustColorscale();

//        refresh = true;
//        adjustColorscale();
//        repaintModel();
    }

	private void checkSetSliderValues(float minValue, float maxValue)
	{
		if(minValue == minSliderValue && maxValue == maxSliderValue) return;
		scale = (maxValue - minValue)/maxIndex;
		minSliderValue = minValue;
		maxSliderValue = maxValue;
		if(debug) System.out.println("ColorscalePanel.checkSetSliderValues:  minSliderValue set to " + minSliderValue + " maxSliderValue set to " + maxSliderValue);
		if(colorscale == null) return;
		if(!colorscale.setEditRange(minSliderValue, maxSliderValue)) return;
//		adjustMinSliderValue(minSliderValue);
//		adjustMaxSliderValue(maxSliderValue);
		adjustColorscale();
	}

    public void setHistogram(float[] data)
    {
        if(data != null)
        {
            histogramPanel.updateData(data, minValue, maxValue);
            hasHistogram = true;
        }
    }

    private void initializeSliderScale()
    {
        if(colorscale != null)
		{
			float[] range = colorscale.getRange();
			minValue = range[0];
			maxValue = range[1];
			float[] editRange = colorscale.getEditRange();
			minSliderValue = editRange[0];
			maxSliderValue = editRange[1];
		}
		if(debug) System.out.println("ColorscalePanel.initializeSliderScale:  minSliderValue set to " + minSliderValue + " maxSliderValue set to " + maxSliderValue);
		maxIndex = nPoints-1;
		scale = (maxValue - minValue)/maxIndex;

	    minSliderIndex = getIndexFromValue(minSliderValue);
	    maxSliderIndex = getIndexFromValue(maxSliderValue);
	    checkSetSliderValues(minSliderValue, maxSliderValue);
        maxSliderBean.setRange(minValue, maxValue);
	    minSliderBean.setRange(minValue, maxValue);
//        maxSliderBean.classInitialize(this, "maxSliderValue", minValue, maxValue, "Max");
//	    minSliderBean.classInitialize(this, "minSliderValue", minValue, maxValue, "Min");


    }

	public int getIndexFromValue(float value)
	{
		int index = Math.round((value - minValue)/scale);
		return StsMath.minMax(index, 0, maxIndex);
	}

    private void initializeSliderPolygons()
    {
        int totalSliderHeight = sliderHeight+heightPerInterval;
        p = new Point[]
        {
            new Point(-sliderWidth, 0),
            new Point(-sliderWidth -(sliderHeight-sliderWidth/2), totalSliderHeight),
            new Point((sliderHeight-sliderWidth/2), totalSliderHeight),
            new Point(0, 0)
        };

        if(!leftSide) flipX = -1;

        minSliderPolygon = new Polygon();
        minSliderPolygon.addPoint(0, 0);
        for(int n = 0; n < 4; n++)
            minSliderPolygon.addPoint(flipX*p[n].x, p[n].y);

        maxSliderPolygon = new Polygon();
        maxSliderPolygon.addPoint(0, heightPerInterval);
        for(int n = 0; n < 4; n++)
            maxSliderPolygon.addPoint(flipX*p[n].x, -p[n].y+heightPerInterval);

        int totalKeyHeight = keyHeight+heightPerInterval;
        p = new Point[]
        {
            new Point((keyHeight-keyWidth/2), -totalKeyHeight/2),
            new Point(-keyWidth, -totalKeyHeight/2),
            new Point((keyHeight-keyWidth/2)-keyWidth, 0),
            new Point(-keyWidth, totalKeyHeight/2),
            new Point((keyHeight-keyWidth/2), totalKeyHeight/2),
            new Point(0, 0)
        };

        keyPolygon = new Polygon();
        keyPolygon.addPoint(0, 0);
        for (int n = 0; n < 6; n++)
            keyPolygon.addPoint(-p[n].x, p[n].y);
    }

    private void drawVoxelKeys(Graphics g)
    {
        int minIdx, maxIdx, x, y, x1, y1, y2;
		StsColorscale.VoxelKey[] voxelKeys = colorscale.getVoxelKeys();
        if((voxelKeys == null) || (!voxelBtn.isSelected())) return;

        x = cbarRect.x;
        y = yOrigin + nPoints;
        g.translate(x, y);
        for(int n = 0; n < voxelKeys.length; n++)
        {
            StsColorscale.VoxelKey key = voxelKeys[n];
            Polygon voxelPolygon = new Polygon();
            if(key != null)
            {
                g.setColor(Color.BLUE);
                x1 = -5;
                minIdx = key.getMin();
                y1 = - minIdx*heightPerInterval;
                voxelPolygon.addPoint(x1, y1);
                maxIdx = key.getMax();
                y2 = - maxIdx*heightPerInterval;
                voxelPolygon.addPoint(x1, y2);
                voxelPolygon.addPoint(x1+10, y2);
                voxelPolygon.addPoint(x1+10, y1);
                voxelPolygon.addPoint(x1, y1);
                g.fillPolygon(voxelPolygon);
                g.setColor(SystemColor.controlDkShadow);
                drawHighlight(g, x1, y1, x1+10, y1, 0, -1);
                drawHighlight(g, x1+10, y1, x1+10, y2, 1, 0);
                g.setColor(SystemColor.controlLtHighlight);
                drawHighlight(g, x1, y1, x1, y2, 1, 0);
                drawHighlight(g, x1, y2, x1+10, y2, 0, -1);
            }
        }
        g.translate(-x, -y);
    }

    private void drawKey(Graphics g) {
        int x, y, i;
        String desc = new String();

        int[] px;
        int[] py;

        g.setColor(Color.black);
        x = cbarRect.x;
        y = yOrigin + nPoints - keyIndex*heightPerInterval;
        g.drawLine(x, y, x+cbarRect.width+4, y);

        px = keyPolygon.xpoints;
        py = keyPolygon.ypoints;

        desc = "[" + keyIndex + "]";
        g.drawString(desc, x+keyWidth, y+10+keyHeight);

        int start = cbarRect.x + sliderWidth/2;
        g.translate(start, y);

        g.setColor(key);
        g.fillPolygon(keyPolygon);

        g.setColor(SystemColor.controlDkShadow);
        drawHighlight(g, px[4], py[4], px[5], py[5], 0, -1);
        drawHighlight(g, px[3], py[3], px[4], py[4], 1, -1);
        g.drawLine(px[5], py[5], px[6], py[6]);
        g.drawLine(px[2], py[2], px[3], py[3]);

        g.setColor(SystemColor.controlLtHighlight);
        drawHighlight(g, px[0], py[0], px[1], py[1], -flipX, 1);
        drawHighlight(g, px[1], py[1], px[2], py[2], 0, 1);

        g.translate(-start, -y);
    }

    private void addKeyPopup(int x, int y)
    {
        JPopupMenu tp = new JPopupMenu("Key Popup");
        this.add(tp);
        tp.add(colorChangeKey);
        tp.add(colorResetKey);
        tp.add(increaseWidthKey);
        tp.add(decreaseWidthKey);
        tp.show(this, x, y);
    }

    private void addVoxelKeyPopup(int x, int y)
    {
        JPopupMenu tp = new JPopupMenu("Voxel Key Popup");
        this.add(tp);
        tp.add(vd);
        tp.show(this, x, y);
    }

    private void drawSlider(Graphics g, Polygon sliderPolygon, int sliderIndex, boolean flipVertical)
    {
        int[] px = sliderPolygon.xpoints;
        int[] py = sliderPolygon.ypoints;

        int x = xOrigin;
        int y = yOrigin + nPoints - sliderIndex*heightPerInterval;

        // draw black line from center of bar to slider tip
        g.setColor(Color.black);
        g.drawLine(x, y, x+flipX*cbarRect.width/2, y);
        g.drawLine(hbarRect.x,y,hbarRect.x+flipX*hbarRect.width/2,y);

        // Write text on amount of data clipped
        g.setColor(Color.blue);
        String desc = colorFormat.formatValue(histogramPanel.getBottomPercentageClipped()) + "% clipped";
        g.drawString(desc, hbarRect.x, yOrigin + nPoints + 10);
        desc = colorFormat.formatValue(histogramPanel.getTopPercentageClipped()) + "% clipped";
        g.drawString(desc, hbarRect.x, yOrigin);

        int start = cbarRect.x + sliderWidth/2;
        g.setColor(Color.black);
        g.translate(start, y);
        g.setColor(SystemColor.controlHighlight);
        g.fillPolygon(sliderPolygon);

        // highlight above and left, shadow below
        if(flipVertical)
        {
            g.setColor(SystemColor.controlLtHighlight);
            drawHighlight(g, px[2], py[2], px[3], py[3], 0, -1);
            drawHighlight(g, px[3], py[3], px[4], py[4], -flipX, -1);

            g.setColor(SystemColor.controlDkShadow);
            drawHighlight(g, px[0], py[0], px[1], py[1], 0, 1);
            g.drawLine(px[1], py[1], px[2], py[2]);
        }
        else
        {
            g.setColor(SystemColor.controlLtHighlight);
            drawHighlight(g, px[0], py[0], px[1], py[1], 0, -1);

            g.setColor(SystemColor.controlDkShadow);
            drawHighlight(g, px[1], py[1], px[2], py[2], 1, 0);
            drawHighlight(g, px[2], py[2], px[3], py[3], 0, 1);
            g.drawLine(px[3], py[3], px[4], py[4]);
        }

        g.translate(-start, -y);
    }

    // draw two lines, one offset by dx,dy from the other
    private void drawHighlight(Graphics g, int x0, int y0, int x1, int y1, int dx, int dy)
    {
        g.drawLine(x0, y0, x1, y1);
        g.drawLine(x0+dx, y0+dy, x1+dx, y1+dy);
    }

    private void addMouseListener()
    {
        MouseListener mouseListener = new MouseAdapter()
        {
            public void mousePressed(MouseEvent e)
            {
                int x, y;

                pickedX = e.getX();
                pickedY = e.getY();

                int start = cbarRect.x + sliderWidth/2;
                x = pickedX - start;
                y = pickedY - (yOrigin + nPoints) + 2;
                if(viewOnly) return;

				sliderSelected = NONE;
				int mods = e.getModifiers();
				if(minSliderPolygon.contains(x, y + minSliderIndex*heightPerInterval))
                {
                    sliderSelected = MIN_SLIDER;
                    pickedIndex = minSliderIndex;
                }
                else if(maxSliderPolygon.contains(x, y + maxSliderIndex*heightPerInterval))
                {
                    sliderSelected = MAX_SLIDER;
                    pickedIndex = maxSliderIndex;
                }
                else if (keyPolygon.contains(x, y + keyIndex * heightPerInterval))
				{
                        sliderSelected = KEY_SLIDER;
                        pickedIndex = keyIndex;
						if((e.isShiftDown() && ((mods & InputEvent.BUTTON3_MASK)) != 0) || ((mods & InputEvent.BUTTON2_MASK) != 0))
							addKeyPopup(pickedX, pickedY);
                }
				else if(voxelBtn.isSelected())
				{
					int idx = -y/heightPerInterval;
                    if((idx >= maxSliderIndex) || (idx <= minSliderIndex))
                        return;
					if(isVoxelKeySelected(idx))
					{
						if((e.isShiftDown() && ((mods & InputEvent.BUTTON3_MASK)) != 0) || ((mods & InputEvent.BUTTON2_MASK) != 0))
							addVoxelKeyPopup(pickedX, pickedY);
					}
					else
					{
						System.out.println("Mouse pressed for voxel key at index " + idx);
						float value = colorscale.getDataValue(idx);
						voxelMinBean.setValue(value);
						voxelMaxBean.setValue(value);
					}
                }
            }

			public void mouseReleased(MouseEvent e)
			{
				mouseMotionEvent(e);
//				adjustColorscale();
				sliderSelected = NONE;
				dataValue = nullValue;
				repaint();
//				ActionEvent ae = new ActionEvent(this, e.getID(), "colorscalePanel mouse released", 0, 0);
//				fireActionPerformed(e);
			}
 /*
 			public void mouseExited(MouseEvent e)
			{
				mouseMotionEvent(e);
//				adjustColorscale();
				sliderSelected = NONE;
				dataValue = nullValue;
				ActionEvent ae = new ActionEvent(this, e.getID(), "colorscalePanel mouse released", 0, 0);
				fireActionPerformed(e);
			}
*/
        };
		this.addMouseListener(mouseListener);
    }

    private boolean isVoxelKeySelected(int idx)
    {
		StsColorscale.VoxelKey[] voxelKeys = colorscale.getVoxelKeys();
		if(voxelKeys == null) return false;

        for(int n = 0; n < voxelKeys.length; n++)
        {
            StsColorscale.VoxelKey key = voxelKeys[n];
			int min = key.getMin();
			int max = key.getMax();
            if(idx >= min - 5 && idx <= max + 5)
            {
                sliderSelected = VOXEL_SLIDER;
                currentVoxelKey = key;
				setVoxelRange();
                return true;
            }
        }
        return false;
    }

    private void addMouseMotionListener()
    {
        MouseMotionListener motionListener = new MouseMotionAdapter()
        {
            public void mouseDragged(MouseEvent e)
            {
                mouseMotionEvent(e);
            }
            public void mouseReleased(MouseEvent e)
            {
                 mouseMotionEvent(e);
            }
              public void mouseExited(MouseEvent e)
            {
                 mouseMotionEvent(e);
            }
        };
        this.addMouseMotionListener(motionListener);
    }

    /** called only if mouseDragged or mouseRelease event occurs */
    private void mouseMotionEvent(MouseEvent e)
    {
		int eventID = e.getID();
//		if(eventID == MouseEvent.MOUSE_PRESSED) return;

//		boolean refresh = true;
        int dy = pickedY - e.getY();
        int index = pickedIndex + dy/heightPerInterval;
		index = StsMath.minMax(index, 0, maxIndex);
        if(sliderSelected == NONE)
        {
            pickedY = e.getY();
            int y = pickedY - (yOrigin + nPoints) + 2;
            index = -y/heightPerInterval;
        }
        if(sliderSelected == MIN_SLIDER)
        {
			if(minSliderIndex == index)
            {
                if(eventID == MouseEvent.MOUSE_RELEASED)
                    adjustColorscale();
                return;
            }
            minSliderIndex = index;
			minSliderValue = computeValueFromIndex(minSliderIndex);
			minSliderBean.setValue(minSliderValue);
			if(!adjustMinSliderValue(minSliderValue)) return;
//            if(syncSliders) syncMax();
			adjustColorscale(eventID);
//            refresh = true;
        }
        else if(sliderSelected == MAX_SLIDER)
        {
			if(maxSliderIndex == index)
            {
                if(eventID == MouseEvent.MOUSE_RELEASED)
                    adjustColorscale();
                return;
            }
            maxSliderIndex = index;
            maxSliderValue = computeValueFromIndex(maxSliderIndex);
			maxSliderBean.setValue(maxSliderValue);
			if(!adjustMaxSliderValue(maxSliderValue)) return;
 //           if(syncSliders) syncMin();
			adjustColorscale(eventID);
//           refresh = true;
        }
        else if(sliderSelected == KEY_SLIDER)
        {
            computeNewColors(index);
			adjustColorscale(eventID);
 //           refresh = true;
        }
        else if(sliderSelected == VOXEL_SLIDER)
        {
            pickedY = e.getY();
            int y = pickedY - (yOrigin + nPoints) + 2;
            index = -y/heightPerInterval;

            if(currentVoxelKey == null) return;
            int maxIdx = currentVoxelKey.getMax();
            int minIdx = currentVoxelKey.getMin();
            int width = maxIdx - minIdx;
            int ctrIdx = minIdx + width/2;
            int diff = index - ctrIdx;
            if((maxIdx + diff) > maxIndex)
            {
                maxIdx = maxIndex;
                minIdx = maxIndex - width;
            }
            else if((minIdx + diff) < 0)
            {
                maxIdx = width;
                minIdx = 0;
            }
            else
            {
                maxIdx = maxIdx + diff;
                minIdx = minIdx + diff;
            }
            colorscale.editVoxelKey(currentVoxelKey, minIdx, maxIdx);
            setVoxelRange();
            setVoxelTransparency();
			adjustColorscale(eventID);
//            refresh = true;
        }
        else
        {
            if(pickedX > cbarRect.x + cbarRect.width)
                return;
            // Output the data value associated with the mouse position. SAJ
            if((index < 0) || (colorscale == null) || (index > maxIndex))
            {
                dataValue = nullValue;
                return;
            }
            /*
            else
            {
            	computeDataValue(e, index);
                repaint();
            }
            */
//            refresh = false;
            if(voxelBtn.isSelected())
            {
                int y = pickedY - (yOrigin + nPoints) + 2;
                int idx = -y/heightPerInterval;
                if(idx >= maxSliderIndex)
                    idx = maxSliderIndex;
                else if(idx < minSliderIndex)
                    idx = 0;

                float value = colorscale.getDataValue(idx);
                float currentMinValue = voxelMinBean.getFloatValue();
                float currentMaxValue = voxelMaxBean.getFloatValue();

                voxelStart = -(computeIndexFromValue(currentMinValue)*heightPerInterval) - 2 + (yOrigin + nPoints);
                voxelEnd = -(computeIndexFromValue(currentMaxValue)*heightPerInterval) - 2 + (yOrigin + nPoints);

                //System.out.println("Dragged to new value " + value + " current range " + currentMinValue + " " + currentMaxValue);
                if(value < currentMinValue)
                    voxelMinBean.setValue(value);
                else if(value > currentMaxValue)
                    voxelMaxBean.setValue(value);
            }
        }
    	computeDataValue(e, index);
    	repaint();
	/*
		if(refresh)
		{
			adjustColorscale();
//			repaint();
		}
	*/
   }

   private void computeDataValue(MouseEvent e, int index)
   {
	   offScreenGraphics.dispose();
	   dataValue = colorscale.getDataValue(index);

	   int width = cbarRect.width + cbarRect.x;
	   dataX = e.getX();
	   if((dataWidth + e.getX()) > width)
		   dataX = width - dataWidth;
	   else if(e.getX() < 0)
		   dataX = 1;

	   dataY = e.getY() - fm.getHeight();
	   if((e.getY() - fm.getHeight())  < cbarRect.y)
		   dataY = cbarRect.y + 3;
	   else if(e.getY() > cbarRect.y + cbarRect.height)
		   dataY = cbarRect.y + cbarRect.height - fm.getHeight() - 3;
   }

   private float computeValueFromIndex(int sliderIndex)
   {
	   return sliderIndex * scale + minValue;
   }

   private int computeIndexFromValue(float value)
   {
	   int index = Math.round((value - minValue)/scale);
	   return StsMath.minMax(index, 0, maxIndex);
   }

    private void computeNewColors(int newIndex)
    {
        if(numKeyIndices == 0) numKeyIndices = 1;

        // Change last Key index back
        if(numKeyIndices > 0)
        {
            newColors = colorscale.getOriginalColors();

                // Change current key index to Key Color
            keyIndex = newIndex;
            if ( (keyIndex + numKeyIndices / 2) > 254)
                keyIndex = 254 - numKeyIndices / 2;
            if ( (keyIndex - numKeyIndices / 2) < 0)
                keyIndex = numKeyIndices / 2;

            for (int i = keyIndex - numKeyIndices / 2; i <= keyIndex + numKeyIndices / 2; i++)
                newColors[i] = key;

//            refresh = true;
//            adjustColorscale();
//            repaintModel();
        }
        return;
    }

    /** Panel contains a controlPanel and a sliderPanel. The control panel has colorscale edit buttons
     *  and the sliderPanel contains the colorbarPanel, histogramPanel, and min max text boxes.
     *  If this is viewOnly, there is no control panel.
     *  If the orientation is VERTICAL, the controlPanel is below the sliderPanel.
     *  If HORIZONTAL, the controlPanel is to the left of the sliderPanel.
     */
    private void jbInit() throws Exception
    {
        this.setLayout(new GridBagLayout());

        colorbarPanel.setFont(null);
        insets = colorbarPanel.getInsets();
        Dimension colorbarSize = new Dimension(barWidth+4, barHeight+10);
        colorbarPanel.setMinimumSize(colorbarSize);
        colorbarPanel.setPreferredSize(colorbarSize);
        Dimension histogramSize = histogramPanel.getSize();
        histogramPanel.setMinimumSize(colorbarSize);
        histogramPanel.setPreferredSize(colorbarSize);
/*
        Dimension histogramSize = new Dimension(histogramPanel.getVerticalBarWidth(), barHeight);
        histogramPanel.setMinimumSize(histogramSize);
        histogramPanel.setPreferredSize(histogramSize);
*/
        Dimension size;

        if(viewOnly)
        {
            size = new Dimension(colorbarSize.width + histogramSize.width + 100, colorbarSize.height);

            this.add(colorbarPanel,    new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
                ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 50, 0, 2), 0, 0));
            this.add(histogramPanel,    new GridBagConstraints(1, 0, 1, 1, 1.0, 1.0
                ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 2, 0, 2), 0, 0));
        }
        else
        {
//            compressBtn.setSelected(true);
            compressBtn.setText("Compress");
            compressBtn.setFont(defaultFont);
            compressBtn.addActionListener(this);

            clipBtn.setText("Clip");
            clipBtn.setFont(defaultFont);
            clipBtn.addActionListener(this);

            flipBtn.setFont(defaultFont);
            flipBtn.addActionListener(this);

            syncBtn.setFont(defaultFont);
            syncBtn.addActionListener(this);

            voxelBtn.setFont(defaultFont);
            voxelBtn.addActionListener(this);

//            transparencyBtn.setSelected(true);
            transparencyBtn.setText("Transparency");
            transparencyBtn.setFont(defaultFont);
            transparencyBtn.addActionListener(this);

            opaqueBtn.setText("Opaque");
            opaqueBtn.setFont(defaultFont);
            opaqueBtn.addActionListener(this);
/*
            maxData.setBorder(BorderFactory.createEtchedBorder());
            maxData.setColumns(10);
            maxData.setHorizontalAlignment(minData.RIGHT);
            maxData.addFocusListener(this);
            maxData.addActionListener(this);

            minData.setBorder(BorderFactory.createEtchedBorder());
            minData.setColumns(10);
            minData.setHorizontalAlignment(minData.RIGHT);
            minData.addFocusListener(this);
            minData.addActionListener(this);

            voxelMaxText.setColumns(5);
			voxelMaxText.setSelectionStart(0);
            voxelMaxText.setHorizontalAlignment(JTextField.RIGHT);
            voxelMinText.setColumns(5);
			voxelMinText.setSelectionStart(0);
            voxelMinText.setHorizontalAlignment(JTextField.RIGHT);
 */
            voxelSetBtn.setText("Set");
            voxelSetBtn.setFont(defaultFont);
            voxelSetBtn.addActionListener(this);

            JPanel sliderPanel = new JPanel();
            //sliderPanel.setLayout(new GridBagLayout());

            JPanel controlPanel = new JPanel();
            controlPanel.setLayout(new GridBagLayout());
            controlPanel.setBorder(BorderFactory.createEtchedBorder());

            JPanel voxelPanel = new JPanel();
            voxelPanel.setLayout(new GridBagLayout());
            voxelPanel.setBorder(BorderFactory.createEtchedBorder());

            if(orientation == StsParameters.VERTICAL)
            {
                size = new Dimension(colorbarSize.width + histogramSize.width + 150, colorbarSize.height + 2 * sliderHeight + 65 + (6 * 20));

                this.add(sliderPanel, new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
                    , GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(2, 2, 2, 2), 0, 0));
                this.add(controlPanel, new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
                    , GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(2, 2, 2, 2), 0, 0));
                this.add(voxelPanel, new GridBagConstraints(0, 2, 1, 1, 1.0, 1.0
                    , GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(2, 2, 2, 2), 0, 0));
            }
            else
            {
                size = new Dimension(colorbarSize.width + histogramSize.width + 300, colorbarSize.height + 2*sliderHeight);

				this.add(sliderPanel, new GridBagConstraints(0, 0, 1, 2, 1.0, 1.0
					, GridBagConstraints.WEST, GridBagConstraints.BOTH,new Insets(2, 2, 2, 2),0, 0));
                this.add(controlPanel, new GridBagConstraints(1, 0, 1, 1, 1.0, 1.0
                    , GridBagConstraints.EAST, GridBagConstraints.BOTH, new Insets(2, 2, 2, 2), 0, 0));
                this.add(voxelPanel, new GridBagConstraints(1, 1, 1, 1, 1.0, 1.0
                    , GridBagConstraints.EAST, GridBagConstraints.BOTH, new Insets(2, 2, 2, 2), 0, 0));
            }

 //           sliderPanel.add(maxLabel, new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
 //               , GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(2, 2, 2, 2), 0, 0));
            //sliderPanel.add(maxSliderBean, new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            //    , GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(2, 2, 2, 2), 0, 0));

            if(!leftSide)
            {
                sliderPanel.add(colorbarPanel, new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0,
                    GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 50, 0, 0), 0, 0));
                sliderPanel.add(histogramPanel, new GridBagConstraints(1, 0, 1, 1, 1.0, 1.0,
                        GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 300, 0, 0), 0, 0));
            }
            else
            {
                sliderPanel.add(colorbarPanel, new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0,
                    GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 50, 0, 0), 0, 0));
                sliderPanel.add(histogramPanel, new GridBagConstraints(1, 1, 1, 1, 1.0, 1.0,
                        GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 2, 0, 2), 0, 0));
            }
 //           sliderPanel.add(minLabel, new GridBagConstraints(0, 2, 1, 1, 1.0, 1.0
 //               , GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(2, 2, 2, 2), 0, 0));
            //sliderPanel.add(minSliderBean, new GridBagConstraints(0, 2, 1, 1, 1.0, 1.0
            //    , GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(2, 2, 2, 2), 0, 0));

            controlPanel.add(compressBtn, new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
                , GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 2, 2, 2), 0, 0));
            controlPanel.add(clipBtn, new GridBagConstraints(1, 0, 1, 1, 1.0, 1.0
                , GridBagConstraints.WEST, GridBagConstraints.NONE,  new Insets(0, 2, 2, 2), 0, 0));
            controlPanel.add(flipBtn, new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
                , GridBagConstraints.WEST, GridBagConstraints.NONE,  new Insets(0, 2, 2, 2), 0, 0));
            controlPanel.add(syncBtn, new GridBagConstraints(1, 1, 1, 1, 1.0, 1.0
                , GridBagConstraints.WEST, GridBagConstraints.NONE,  new Insets(0, 2, 2, 2), 0, 0));
//            controlPanel.add(rotateSlider, new GridBagConstraints(1, 4, 2, 1, 1.0, 1.0
//                , GridBagConstraints.WEST, GridBagConstraints.NONE,  new Insets(2, 2, 2, 2), 0, 0));
            controlPanel.add(transparencyBtn, new GridBagConstraints(0, 2, 1, 1, 1.0, 1.0
                , GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 2, 2, 2), 0, 0));
            controlPanel.add(opaqueBtn, new GridBagConstraints(1, 2, 1, 1, 1.0, 1.0
                , GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 2, 2, 2), 0, 0));

            voxelPanel.add(voxelBtn, new GridBagConstraints(0, 0, 1, 1, 0.0, 1.0
                , GridBagConstraints.WEST, GridBagConstraints.NONE,  new Insets(0, 2, 2, 2), 0, 0));
            voxelPanel.add(voxelMinBean, new GridBagConstraints(1, 0, 1, 1, 1.0, 1.0
                , GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL,  new Insets(0, 2, 2, 2), 0, 0));
            voxelPanel.add(voxelMaxBean, new GridBagConstraints(2, 0, 1, 1, 1.0, 1.0
                , GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL,  new Insets(0, 2, 2, 2), 0, 0));
            voxelPanel.add(voxelSetBtn, new GridBagConstraints(3, 0, 1, 1, 0.0, 1.0
                , GridBagConstraints.WEST, GridBagConstraints.NONE,  new Insets(0, 2, 2, 2), 0, 0));

            modeGroup.add(clipBtn);
            modeGroup.add(compressBtn);
            transparencyGroup.add(opaqueBtn);
            transparencyGroup.add(transparencyBtn);

            if(showVoxelControls == false)
                remove(voxelPanel);
        }
        setMinimumSize(size);
        setPreferredSize(size);
    }

    private void paintOffScreen(Graphics og)
    {
        og.clearRect(0, 0, dataWidth, dataHeight);
        og.setColor(Color.red);
        og.drawRect(0, 0, dataWidth-1, dataHeight-1);
        og.setColor(Color.black);
        og.setFont(floatFont);
        int x = (int)((float)dataWidth * 0.05f);
        int y = (int)((float)dataHeight - (float)dataHeight*0.10f);
        og.drawString(labelFormat.format(dataValue), x, y);
    }

    public void setColorscale(StsColorscale colorscale)
	{
        this.colorscale = colorscale;
        this.mode = COLORSCALE;
        init();
    }
    public StsColorscale getColorscale() { return colorscale; }
    public float[] getHistogram() { return histogramPanel.getDataHistogram(); }

    public StsSpectrum getSpectrum()
    {
        if(mode == COLORSCALE)
            return colorscale.getSpectrum();
        else if(mode == SPECTRUM)
            return this.spectrum;
        else if(mode == BEAN)
            return (StsSpectrum)fieldBean.getPanelObject();
        else
            return null;
    }

    public void setSpectrum(StsSpectrum spectrum)
    {
        if(mode == COLORSCALE)
        {
            colorscale.setSpectrum(spectrum);
 //           refresh = true;
            adjustColorscale();
//            adjustColorscale();
//            repaintModel();
        }
        else
        {
            this.spectrum = spectrum;
            repaint();
        }
    }

    public boolean setColors(Color[] colors)
    {
        if(mode != COLORS) return false;
        this.inColors = colors;
        Graphics gt = getGraphics();
        paint(gt);
        return true;
    }

    public void setRange(float[] range) { colorscale.setRange(range); }
    public void setRange(float min, float max) { colorscale.setRange(min, max); }
    public void setEditRange(float[] range)
    {
        colorscale.setEditRange(range);
    }

    public void setTransparencyMode(boolean mode)
    {
        colorscale.setTransparencyMode(mode);
        adjustColorscale();
    }
    public void setCompressionMode(int mode)
    {
        colorscale.setCompressionMode(mode);
        adjustColorscale();
    }

    public void paint(Graphics g)
    {
        Color[] colors;
        int nColors;
        float[] range;
        float min;
        float max;

        if(g == null) return;
        super.paint(g);

        if(mode == COLORSCALE && colorscale == null)
            return;
        else if( (mode == SPECTRUM || mode == BEAN) && spectrum == null)
            return;

//		minSliderBean.setValue(minSliderValue);
//		maxSliderBean.setValue(maxSliderValue);

        if(offScreenImage == null)
        {
            dataWidth = fm.stringWidth("99999.99");
            dataWidth += dataWidth*0.10;
            dataHeight = fm.getHeight();
            dataHeight += dataHeight*0.10;
            offScreenImage = createImage(dataWidth, dataHeight);
        }
        offScreenGraphics = offScreenImage.getGraphics();

        hbarRect = histogramPanel.getBounds();

        cbarRect = colorbarPanel.getBounds();
        cInsets = colorbarPanel.getInsets();
        if(mode == SPECTRUM) cbarRect = this.getBounds();

        xOrigin = cbarRect.x;
        if ((!leftSide) && (!viewOnly))
            xOrigin += cbarRect.width;
        yOrigin = hbarRect.y;

        back = this.getVisibleRect();
        g.setColor(getBackground());
        g.fillRect(back.x, back.y, back.width, back.height);
        paintChildren(g);

//        reinitializeSliderScale();
        histogramPanel.setClip(minSliderIndex, maxSliderIndex);

        //SAJ Some Cards will hang, Need to make this occur only when the idx has changed
//        if(hasHistogram) histogramPanel.repaint();

        if(mode == COLORSCALE)
        {
            if(colorscale == null) return;
            colors = colorscale.getNewColors();
            nColors = colorscale.getNColors();
            range = colorscale.getRange();
            min = range[0];
            max = range[1];
        }
        else if(mode == SPECTRUM)
        {
            if(spectrum == null) return;
            colors = spectrum.getColors();
            nColors = spectrum.getNColors();
            min = 0.0f;
            max = 100.0f;
        }
        else if(mode == BEAN)
        {
            if(spectrum == null) return;
            spectrum = (StsSpectrum)fieldBean.getPanelObject();
            colors = spectrum.getColors();
            nColors = spectrum.getNColors();
            min = 0.0f;
            max = 100.0f;
        }
        else
        {
            colors = inColors;
            nColors = colors.length;
            min = 0.0f;
            max = 100.0f;
        }

        if(nColors > 0)
		{
            int maxLabelWidth = 0;
            // determine how big each box should be
            g.setFont(new Font("Dialog",0,10));
            FontMetrics fm = g.getFontMetrics();
            int fontSize = fm.getAscent() + fm.getHeight();

            int height = heightPerInterval * (255/nColors);
            int totalHeight = nColors;
            int width = back.width;

            if(mode == COLORSCALE)
            {
                // draw range labels if there is enough room
                int numLabels = nColors;
                int maxLabels = (height*255) / fontSize;

                if (maxLabels < numLabels) {
                    if (maxLabels < 2)
                        numLabels = 0;
                    else if (maxLabels > 11)
                        numLabels = 11;
                    else if (maxLabels > 9)
                        numLabels = 9;
                    else if (maxLabels > 7)
                        numLabels = 7;
                    else if (maxLabels > 5)
                        numLabels = 5;
                    else if (maxLabels > 3)
                        numLabels = 3;
                    else
                        numLabels = 2;
                }
                int rectWidth = width;

                // Draw Labels
                if ((numLabels > 0) && (labelsOn))
                {
                    float labelSpacing = (float)(height * 255) / (float)(numLabels - 1);
                    float inc = (max - min) / (numLabels - 1);
                    colorFormat.formatValueRange(min, max);
                    String[] values = new String[numLabels];
                    for (int i = 0; i < numLabels; i++)
                    {
                        float value = min + (inc * i);
                        values[i] = colorFormat.format(value);
                        int labelWidth = fm.stringWidth(values[i]);
                        if (labelWidth > maxLabelWidth)
                            maxLabelWidth = labelWidth;
                    }
                    if (maxLabelWidth < rectWidth)
                        rectWidth = width - maxLabelWidth - 10;

                    int xpos = cbarRect.x -10 - maxLabelWidth;
                    if(xpos < 0) xpos = 0;

                    for (int i = 0; i < numLabels; i++)
                    {
                        g.setColor(Color.black);
                        int ypos = yOrigin + (int)(totalHeight - (i * labelSpacing));
                        g.drawString(values[i], xpos, ypos);
                        g.drawLine(cbarRect.x, ypos, cbarRect.x-5, ypos);
                    }
                }
            }
            // Draw the Colors
            int x = xOrigin + 1;
            if ((!leftSide) && (!viewOnly))
                x = xOrigin - cbarRect.width + 1;
            int y = yOrigin;

            if(mode != COLORSCALE) x = 1;

            // Draw Color Scale - Bottom up

			//System.out.println("Panel " + getName() + " Drawing colors: botColor (0) " + colorToString(colors[0]) + " topColor ("  + nColors + "-1)" + colorToString(colors[nColors-1]));
            //TODO bug here: length of colors is 255 while nColors is 256
            for (int i = nColors - 1; i >= 0; i--)
            {
                // draw a rectangle for each color
                g.setColor(colors[nColors - 1 - i]);
                g.fillRect(x, y + i*height, cbarRect.width + 1, height);
            }
            g.setColor(Color.black);
            g.drawRect(x-1, y, cbarRect.width + 2, heightPerInterval*nColors);
        }

        if (!viewOnly) {
            // Draw Slider Keys
            drawSlider(g, minSliderPolygon, minSliderIndex, false);
            drawSlider(g, maxSliderPolygon, maxSliderIndex, true);

            // Draw Dynamic Contour Key
            drawKey(g);

            // Draw Voxel Keys
            drawVoxelKeys(g);
            if(voxelBtn.isSelected())
            {
                g.setColor(Color.BLACK);
                g.drawLine(cbarRect.x - 10,voxelStart,cbarRect.x + 10,voxelStart);
                g.drawLine(cbarRect.x - 10,voxelEnd,cbarRect.x + 10,voxelEnd);
            }
        }
//        if(mode == COLORSCALE)
//            adjustColorscale();

        if(dataValue != StsParameters.nullValue)
        {
            paintOffScreen(offScreenGraphics);
            g.drawImage(offScreenImage, dataX, dataY, this);
        }
    }

	private String colorToString(Color color)
	{
	 return getClass().getName() + "[r=" + color.getRed() + ",g=" + color.getGreen() + ",b=" + color.getBlue()  + ",a=" + color.getAlpha() + "]";

	}
/*
    public void adjustSpectrum(boolean refresh)
    {
        refresh = true;
        adjustSpectrum();
        repaint();
    }
*/
     private void initializeColorscale()
    {
        if(colorscale == null) return;
		colorscale.adjustSpectrum();
//		if(!colorscale.setEditRange(minSliderValue, maxSliderValue)) return;
        if(numKeyIndices > 0)
            colorscale.setSomeNewColors(keyIndex - numKeyIndices/2, numKeyIndices, key);
		repaintModel();
		colorscale.resetOpacity();
    }

    private void adjustColorscale(int mouseMotionID)
    {
        if(colorscale == null) return;
		colorscale.adjustSpectrum();
//		if(!colorscale.setEditRange(minSliderValue, maxSliderValue)) return;
        if(numKeyIndices > 0)
            colorscale.setSomeNewColors(keyIndex - numKeyIndices/2, numKeyIndices, key);
		repaintModel();
		colorscale.resetOpacity();
		colorscale.colorsChanged(mouseMotionID);
    }

    private void adjustColorscale()
    {
        if(colorscale == null) return;
		colorscale.adjustSpectrum();
//		if(!colorscale.setEditRange(minSliderValue, maxSliderValue)) return;
        if(numKeyIndices > 0)
            colorscale.setSomeNewColors(keyIndex - numKeyIndices/2, numKeyIndices, key);
		repaintModel();
		colorscale.resetOpacity();
		colorscale.colorsChanged(MouseEvent.MOUSE_RELEASED);
    }

    public void repaintModel()
    {
		paint(getGraphics());
		if(observerColorscalePanel != null) observerColorscalePanel.repaint();
    }

    public void setViewOnly(boolean tf)
    {
        viewOnly = tf;
    }

    public void setLabelsOn(boolean lo)
    {
        labelsOn = lo;
    }

    public Dimension getPreferredSize()
    {
        int nColors;
        Dimension dim = super.getPreferredSize();
        if(colorscale != null)
            nColors = colorscale.getNColors();
        else if(spectrum != null)
            nColors = spectrum.getNColors();
        else
            nColors = inColors.length;
        if( dim.height < nColors * perColorHeight )
        {
            FontMetrics fm = getGraphics().getFontMetrics();
            int fontSize = fm.getAscent() + fm.getHeight();
            dim.height = nColors * perColorHeight + fontSize;
        }
        return dim;

    }
/*
    private void computeMinSliderValue()
    {
        minSliderValue = minSliderIndex*scale + minValue;
    }

    private void computeMaxSliderValue()
    {
        maxSliderValue = maxSliderIndex*scale + minValue;
    }
*/

    public void setMaxSliderValue(float value)
	{
		if(colorscale.getEditMax() == value) return;
//        maxSliderValue = value;
		adjustMaxSliderValue(value);
	}

	private boolean adjustMaxSliderValue(float value)
	{
//        float[] range = colorscale.range();
        maxSliderIndex = computeIndexFromValue(value);
		maxSliderValue = computeValueFromIndex(maxSliderIndex);
		if(debug) System.out.println("ColorscalePanel.adjustMaxSliderValue:  maxSliderValue set to " + maxSliderValue);
		if(syncSliders) syncMin();
		if(colorscale == null) return false;
		return colorscale.setEditMax(maxSliderValue);
//		refresh = true;
//		adjustColorscale();
//		paint(getGraphics());
//		repaintModel();
    }

	public float getMaxSliderValue()
	{
        if(colorscale == null) return 0.0f;
		maxSliderValue = colorscale.getEditMax();
		return maxSliderValue;
	}

	public float getMinSliderValue()
	{
        if(colorscale == null) return 0.0f;
		minSliderValue = colorscale.getEditMin();
		return minSliderValue;
	}

    private void syncMin()
    {
        minSliderIndex = maxIndex - maxSliderIndex;
        if(minSliderIndex >= maxSliderIndex)
            minSliderIndex = maxSliderIndex-1;
        else if(minSliderIndex < 0)
            minSliderIndex = 0;
		minSliderValue = minSliderIndex*scale + minValue;
		if(debug) System.out.println("ColorscalePanel.syncMin:  minSliderValue set to " + minSliderValue);
		minSliderBean.setValue(minSliderValue);
		colorscale.setEditMin(minSliderValue);
//		computeMinSliderValue();
    }

	private boolean sliderChanged(float newMin, float newMax)
	{
		return newMin != minSliderValue || newMax != maxSliderValue;
	}

	public void setMinSliderValue(float value)
    {
		if(minSliderValue == value) return;
//        minSliderValue = value;
		adjustMinSliderValue(value);
    }

	private boolean adjustMinSliderValue(float value)
	{
		minSliderIndex = computeIndexFromValue(value);
		minSliderValue = computeValueFromIndex(minSliderIndex);
		if(debug) System.out.println("ColorscalePanel.adjustMinSliderValue:  minSliderValue set to " + minSliderValue);
		if (syncSliders) syncMax();
		if(colorscale == null) return false;
		return colorscale.setEditMin(minSliderValue);
//		refresh = true;
//		adjustColorscale();
//		paint(getGraphics());
//		repaintModel();
    }

    private void syncMax()
    {
		maxSliderIndex = maxIndex - minSliderIndex;
		if(maxSliderIndex <= minSliderIndex)
			maxSliderIndex = minSliderIndex+1;
		else if(maxSliderIndex > maxIndex)
			maxSliderIndex = maxIndex;
		 maxSliderValue = maxSliderIndex*scale + minValue;
		 colorscale.setEditMax(maxSliderValue);
		 maxSliderBean.setValue(maxSliderValue);
		 if(debug) System.out.println("ColorscalePanel.syncMax:  maxSliderValue set to " + maxSliderValue);
    }

    public int getHeightPerInterval() { return heightPerInterval; }
    public void setHeightPerInterval(int heightPerInterval) { this.heightPerInterval = heightPerInterval; }
    public int getBarWidth() {  return barWidth; }
    public void setBarWidth(int barWidth) { this.barWidth = barWidth; }
    public int getNPoints() { return nPoints; }
    public void setNPoints(int nPoints) { this.nPoints = nPoints; }

    public void stateChanged(ChangeEvent e)
    {
        int i;
        Object source = e.getSource();

        if (source == rotateSlider)
        {
            colorscale.setRotateAmount(rotateSlider.getValue());
        }
    }

    public void actionPerformed(ActionEvent e)
    {
        int i;
        Object source = e.getSource();
        // Process Radio Button Input
        if (source == clipBtn)
        {
            setCompressionMode(colorscale.CLIPPED);
        }
        if (source == voxelBtn)
        {
            if(voxelBtn.isSelected())
            {
                voxelMinBean.setEnabled(true);
                voxelMaxBean.setEnabled(true);
                voxelSetBtn.setEnabled(true);
                setVoxelTransparency();
            }
            else
            {
                voxelMinBean.setEnabled(false);
                voxelMaxBean.setEnabled(false);
                voxelSetBtn.setEnabled(false);
                resetVoxelTransparency();
            }
        }
        else if (source == voxelSetBtn)
        {
            // Record Voxel Settings
            float value = voxelMinBean.getFloatValue();
            int idxMin = getIndexFromValue(value);
            value = voxelMaxBean.getFloatValue();
            int idxMax = getIndexFromValue(value);

			if(currentVoxelKey != null)
			{
				if(currentVoxelKey.min == idxMin && currentVoxelKey.max == idxMax)
					return;
			}
            currentVoxelKey = colorscale.addVoxelKey(idxMin, idxMax);
            StsColorscale.VoxelKey[] voxelKeys = colorscale.getVoxelKeys();
            voxelStart = 0;
            voxelEnd = 0;
            setVoxelTransparency();
        }
        else if (source == vd)
        {
            currentVoxelKey = colorscale.deleteVoxelKey(currentVoxelKey);
			if(colorscale.canDrawVoxels())
				setVoxelTransparency();
			else
				resetVoxelTransparency();
            setVoxelRange();
        }

        else if (source == flipBtn)
        {
            colorscale.setFlip(flipBtn.isSelected());
			repaintModel();
        }
        else if(source == syncBtn)
        {
			boolean synched = syncBtn.isSelected();
			colorscale.setIsSynched(synched);
            syncSliders = synched;
        }
        else if (source == compressBtn)
        {
            setCompressionMode(colorscale.COMPRESSED);
        }
        else if (source == transparencyBtn)
        {
            setTransparencyMode(true);
        }
        else if (source == opaqueBtn)
        {
            setTransparencyMode(false);
        }
	/*
        else if (source == minData)
        {
            setMinSliderFromValue(Float.valueOf(minData.getText().trim()).floatValue());
            refresh = true;
        }
        else if (source == maxData)
        {
            setMaxSliderFromValue(Float.valueOf(maxData.getText().trim()).floatValue());
            refresh = true;
        }
	*/
        else if (source == colorChangeKey)
        {
            if(jColorChoose == null)  jColorChoose = new JColorChooser();
            key = jColorChoose.showDialog(this,"Color Selection", key);
            if(numKeyIndices == 0)
                numKeyIndices = 1;
            for(i=keyIndex - numKeyIndices/2; i<=keyIndex + numKeyIndices/2; i++)
                newColors[i] = key;
            adjustColorscale();
        }
        else if (source == colorResetKey)
        {
            colorscale.setNewColors(colorscale.getOriginalColors());
            newColors = (Color []) StsMath.arraycopy(colorscale.getOriginalColors(), colorscale.getNColors());
            keyIndex = newColors.length/2;
            key = newColors[keyIndex];
            numKeyIndices = 0;
            adjustColorscale();
        }
        else if (source == increaseWidthKey)
        {
            if(numKeyIndices == 0)
            {
                newColors[keyIndex] = key;
                numKeyIndices = 1;
            }
            else
            {
                newColors[keyIndex + numKeyIndices / 2 + 1] = key;
                newColors[keyIndex - numKeyIndices / 2 - 1] = key;
                numKeyIndices = numKeyIndices + 2;
            }
            adjustColorscale();
        }
        else if (source == decreaseWidthKey)
        {
            if(numKeyIndices > 1)
            {
                newColors[keyIndex + numKeyIndices/2] = colorscale.getOriginalColors()[keyIndex + numKeyIndices/2];
                newColors[keyIndex - numKeyIndices/2] = colorscale.getOriginalColors()[keyIndex + numKeyIndices/2];
                numKeyIndices = numKeyIndices - 2;
            }
            adjustColorscale();
        }
    }

    public void setOpacityValues(float[] opacity)
    {
        int nColors = opacity.length + 1;
        int[] colorscaleOpacityValues = new int[nColors];
        for(int n = 0; n < nColors-1; n++)
            colorscaleOpacityValues[n] = (int)(opacity[n]*255);
        colorscale.setOpacityValues(colorscaleOpacityValues);
 //       refresh = true;
        adjustColorscale();
//        repaintModel();
    }

    public void setVoxelTransparency()
    {
		//StsColorscale.VoxelKey[] keys = colorscale.getVoxelKeys();
		//if(keys == null || keys.length == 0) return;
        float opacity[] = new float[255];
        for(int i=0; i< 255; i++)
            opacity[i] = 0.0f;
		StsColorscale.VoxelKey[] keys = colorscale.getVoxelKeys();
		if(keys == null || keys.length == 0) return;
		for(int i=0; i<keys.length; i++)
        {
            for(int j=keys[i].getMin(); j<keys[i].getMax(); j++)
                opacity[j] = 1.0f;
        }
        setOpacityValues(opacity);
        histogramPanel.resetOpacityValues();
        return;
    }

    public void resetVoxelTransparency()
    {
        float opacity[] = new float[255];
        for(int i=0; i< 255; i++)
            opacity[i] = 1.0f;
        setOpacityValues(opacity);
        histogramPanel.resetOpacityValues();
    }
    public void focusLost(FocusEvent e) { }
    public void focusGained(FocusEvent e) { }

    public void addActionListener(ActionListener actionListener)
    {
        actionListeners = (ActionListener[])StsMath.arrayAddElement(actionListeners, actionListener);
    }

    public void deleteActionListener(ActionListener actionListener)
    {
        actionListeners = (ActionListener[])StsMath.arrayDeleteElement(actionListeners, actionListener);
    }

    protected void fireActionPerformed(ActionEvent e)
    {
        if(actionListeners == null) return;
        for(int n = 0; n < actionListeners.length; n++)
            actionListeners[n].actionPerformed(e);
    }

    protected void fireActionPerformed(MouseEvent e)
    {
        if(actionListeners == null || actionListeners.length == 0) return;
        ActionEvent ae = new ActionEvent(this, ActionEvent.ACTION_PERFORMED, "colorscalePanel mouse event", 0, 0);
        for(int n = 0; n < actionListeners.length; n++)
            actionListeners[n].actionPerformed(ae);
    }
	static public void main(String[] args)
	 {
		 try
		 {
			 ObjectPanelTest objectPanelTest = new ObjectPanelTest();
             StsModel model = StsModel.constructor();
             StsObjectPanel panel = StsObjectPanel.constructor(objectPanelTest, true);
			 StsToolkit.createDialog(panel, false);
//            objectPanelTest.printState();
//			panel = new StsObjectPanel(objectPanelTest, true);
//			StsToolkit.createDialog(panel);
//			objectPanelTest.printState();
		 }
		 catch(Exception e) { e.printStackTrace(); }
    }
/*
    static public void main(String[] args)
    {
    	Color[] colors =
        {
            Color.blue,
            Color.cyan,
            Color.green,
            Color.magenta,
            Color.orange,
            Color.pink,
            Color.red,
            Color.yellow
        };

        StsColor[] stsColors =
        {
            new StsColor(Color.blue),
            new StsColor(Color.cyan),
            new StsColor(Color.green),
            new StsColor(Color.magenta),
            new StsColor(Color.orange),
            new StsColor(Color.pink),
            new StsColor(Color.red),
            new StsColor(Color.yellow)
        };

        try
        {
            StsModel model = new StsModel();
//            StsSpectrum spectrum = new StsSpectrum("test", colors);
            StsSpectrum spectrum = StsSpectrum.constructor("test", stsColors, 255);
            StsColorscale colorscale = new StsColorscale(spectrum, -100.0f, 100.0f);
            StsColorscalePanel cs = new StsColorscalePanel(colorscale, false);
			StsToolkit.createDialog(cs, true, 200, 500);
            // now rerun it with changes in colorscale
			StsToolkit.createDialog(cs, true, 200, 500);
        }
        catch(Exception e)
        {
            e.printStackTrace();
        }
    }
*/
}

