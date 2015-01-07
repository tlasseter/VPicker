package com.Sts.Actions.Wizards.Color;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;
import java.text.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsSpectrumDialog extends JDialog implements WindowListener, ActionListener, MouseListener, ChangeListener, FocusListener
{
    private ButtonGroup btnGroup1 = new ButtonGroup();
    private ButtonGroup btnGroup2 = new ButtonGroup();

    public StsModel model = null;
    private boolean saved = true;
    private boolean minSet = false;
    private boolean maxSet = false;

    private Frame frame = new Frame();

    JPanel jPanel1 = new JPanel();
    JPanel mainPanel = new JPanel();
    JPanel colorscalePanel = new JPanel();

    private JMenuItem ka = new JMenuItem("Add Key");
    private JMenuItem kc = new JMenuItem("Change Key Color");
    private JMenuItem kd = new JMenuItem("Delete Key");
    private JMenuItem kr = new JMenuItem("Reset Spectrum");

    JPanel btnPanel = new JPanel();
    JButton resetBtn = new JButton();
    JButton saveAsBtn = new JButton();
    JButton previewBtn = new JButton();
    JButton cancelBtn = new JButton();

    JColorChooser jColorChoose = new JColorChooser();

    JPanel modePanel = new JPanel();
    JRadioButton multipleClr = new JRadioButton();
    JRadioButton keyClr = new JRadioButton();

    JPanel paintModePanel = new JPanel();
    JButton activeClr = new JButton();
    JRadioButton paintClr = new JRadioButton();
    JRadioButton indexClr = new JRadioButton();

    StsColor[] keys = null;
    JRadioButton interpolateClrs = new JRadioButton();

    JPanel colorSelect = new JPanel();

    GridBagLayout gridBagLayout6 = new GridBagLayout();
    GridBagLayout gridBagLayout1 = new GridBagLayout();
    GridBagLayout colorSelectLayout = new GridBagLayout();
    GridBagLayout modeLayout = new GridBagLayout();
    GridBagLayout interpolateLayout = new GridBagLayout();

    Font defaultFont = new Font("Dialog",0,11);
    Font graphicFont = new Font("Dialog",0,9);
    StsCursor cursor = null;

    private static final int SINGLE_SELECT = 0;
    private static final int MULTIPLE_SELECT = 1;
    private static final int KEY_SELECT = 2;
    int mode = SINGLE_SELECT;

    boolean isOK = false;

    Insets insets = null;
    int nPoints = 255;
    int maxIndex = 254;
    int heightPerInterval = 1;
    int barHeight;
    int barWidth = 25;
    int minSliderIndex;
    int maxSliderIndex;
    int startIndex = 0;
    Rectangle barRect;
    int xOrigin, yOrigin;
    Polygon minSliderPolygon, maxSliderPolygon, keyPolygon[];
    int pickedX, pickedY;
    int pickedIndex = 0;
    int activeIndex = 0;
    int flipX = 1;
    float maxValue = 1.0f;
    float minValue = 0.0f;
    float maxSliderValue = 0.25f;
    float minSliderValue = 0.75f;
    float scale;
    boolean success = false;
    float[] range = new float[2];

    byte sliderSelected = NONE;
    byte activeSlider = NONE;
    static private byte NONE = 0;
    static private byte MIN_SLIDER = -1;
    static private byte MAX_SLIDER = 1;

    static int sliderHeight = 10;
    static int sliderWidth = 15;
    static int keyHeight = 6;
    static int keyWidth = 10;

    int selectedKeyIndex = -1;
    int keyIndex = -1;
    boolean paintbrushMode = false;

    Point[] p;
    Color[] newColors = null;
    Color activeColor = new Color(255,0,255,255);
    Color nullColor = new Color(255,255,255,255);
    StsSpectrum spectrum = null;
    StsColorscale colorscale = null;
    StsColorscalePanel origColorscalePanel = null;

    boolean colorscaleAvailable = false;

    DecimalFormat labelFormat = new DecimalFormat("#,##0.0###");
    JPanel jPanel2 = new JPanel();
    JLabel jLabel1 = new JLabel();
    JButton nullClrBtn = new JButton();
    JLabel jLabel2 = new JLabel();
    JLabel jLabel3 = new JLabel();
    JLabel jLabel4 = new JLabel();

    JTextField nullRed = new JTextField();
    JTextField nullGreen = new JTextField();
    JTextField nullBlue = new JTextField();
    JTextField nullAlpha = new JTextField();

    JLabel jLabel5 = new JLabel();
    GridBagLayout gridBagLayout2 = new GridBagLayout();
    GridBagLayout gridBagLayout3 = new GridBagLayout();

    public StsSpectrumDialog()
    {
    }

    public StsSpectrumDialog(StsSpectrum spectrum)
    {
        this();
        this.setSpectrum(spectrum);
        nullColor = spectrum.getTransparentColor().getColor();
    }

    public StsSpectrumDialog(Frame frame, String title, StsModel model, boolean modal)
    {
        super(frame, title, modal);
        this.setLocationRelativeTo(frame);
        if(frame != null)
            this.frame = frame;
        this.model = model;

        this.setTitle(title);
        newColors = new Color[255];
        newColors[0] = new Color(255,255,255,255);
        newColors[254] = new Color(0,0,0,255);
        interpolateClrs(0, 254);
        nullColor = new Color(255,255,255,255);
        keys = new StsColor[3];
        keys[0] = new StsColor(newColors[0]);
        keys[0].idx = 0;
        keys[1] = new StsColor(newColors[128]);
        keys[1].idx = 128;
        keys[2] = new StsColor(newColors[254]);
        keys[2].idx = 254;
        range[0] = 0;
        range[1] = newColors.length-1;
        initialize();

        try
        {
            jbInit();
            pack();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }

    }

    public StsSpectrumDialog(String title, StsModel model, StsSpectrum spectrum, boolean modal)
    {
        this.model = model;
        if(model != null) this.frame = model.win3d;
        this.spectrum = spectrum;
        nullColor = spectrum.getTransparentColor().getColor();
        keys = spectrum.getKeys();
        newColors = (Color []) StsMath.arraycopy(spectrum.getColors(), spectrum.getNColors());
        range[0] = 0;
        range[1] = newColors.length-1;
        initialize();

        try
        {
            jbInit();
            pack();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public StsSpectrumDialog(String title, StsModel model, StsColorscale colorscale, boolean modal)
    {
        this.model = model;
        this.setTitle(title);
        if(model != null) this.frame = model.win3d;
        this.spectrum = colorscale.getSpectrum();
        this.colorscale = colorscale;
        nullColor = spectrum.getTransparentColor().getColor();
        colorscaleAvailable = true;
        range = colorscale.getRange();
        keys = colorscale.getSpectrum().getKeys();
        newColors = (Color []) StsMath.arraycopy(colorscale.getOriginalColors(), colorscale.getNColors());
        initialize();

        try
        {
            jbInit();
            pack();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public StsSpectrumDialog(String title, StsModel model, StsColorscalePanel colorscalePanel, boolean modal)
    {
        this(title, model, colorscalePanel.getColorscale(), modal);
        this.origColorscalePanel = colorscalePanel;
    }


    public StsSpectrumDialog(StsModel model)
    {
        this(model.win3d,"New Spectrum", model, true);
    }

    private void initialize() {
        addMouseListener();
        addMouseMotionListener();

        initializeSliderScale();
        initializeSliderPolygons();
        initializeKeyPolygons();
        barHeight = maxIndex*heightPerInterval;

        this.cursor = new StsCursor(this,Cursor.DEFAULT_CURSOR);
    }

    private void jbInit() throws Exception
    {
        this.setModal(false);
        this.addWindowListener(this);

        colorscalePanel.setFont(null);
        insets = colorscalePanel.getInsets();
        Dimension colorbarSize = new Dimension(barWidth+insets.left+insets.right + 30, barHeight+insets.top+insets.bottom + 100);
        colorscalePanel.setMinimumSize(colorbarSize);
        colorscalePanel.setPreferredSize(colorbarSize);
        colorscalePanel.setBorder(BorderFactory.createEtchedBorder());
//        colorscalePanel.setLayout(colorscaleLayout);

        jPanel1.setLayout(gridBagLayout3);

        btnPanel.setBorder(BorderFactory.createEtchedBorder());
        btnPanel.setOpaque(true);
        btnPanel.setLayout(gridBagLayout1);

        modePanel.setBorder(BorderFactory.createEtchedBorder());
        modePanel.setOpaque(true);
        modePanel.setLayout(modeLayout);

        multipleClr.setText("Colors");
        multipleClr.addActionListener(this);
        multipleClr.setSelected(true);
        multipleClr.setFont(defaultFont);
        keyClr.setText("Color Keys");
        keyClr.addActionListener(this);
        keyClr.setFont(defaultFont);

        paintModePanel.setBorder(BorderFactory.createEtchedBorder());
        paintModePanel.setOpaque(true);
        paintModePanel.setLayout(modeLayout);

        activeClr.setText("     ");
        activeClr.setBackground(activeColor);
        activeClr.addActionListener(this);
        activeClr.setOpaque(true);
        paintClr.setText("Drag");
        paintClr.setFont(defaultFont);
        paintClr.setSelected(true);
        paintClr.addActionListener(this);
        indexClr.setText("Index");
        indexClr.setFont(defaultFont);
        indexClr.addActionListener(this);

        resetBtn.setText("Reset");
        resetBtn.setFont(defaultFont);
        resetBtn.addActionListener(this);
        saveAsBtn.setText("Save As...");
        saveAsBtn.setFont(defaultFont);
        saveAsBtn.addActionListener(this);
        previewBtn.setText("Preview");
        previewBtn.setFont(defaultFont);
        previewBtn.addActionListener(this);
        cancelBtn.setText("Cancel");
        cancelBtn.setFont(defaultFont);
        cancelBtn.addActionListener(this);

        jPanel2.setBorder(BorderFactory.createEtchedBorder());
        jPanel2.setLayout(gridBagLayout2);
        jLabel1.setText("Null Value:");
        nullClrBtn.setBackground(Color.black);
        nullClrBtn.setText("");
        nullClrBtn.addActionListener(this);
        jLabel2.setText("Red");
        jLabel3.setText("Green");
        jLabel4.setText("Blue");
        nullRed.setText("0");
        nullRed.setScrollOffset(255);
        nullRed.addActionListener(this);
        nullRed.addFocusListener(this);
        nullGreen.setText("0");
        nullGreen.addActionListener(this);
        nullGreen.addFocusListener(this);
        nullBlue.setText("0");
        nullBlue.addActionListener(this);
        nullBlue.addFocusListener(this);
        jLabel5.setText("Alpha");
        nullAlpha.setText("255");
        nullAlpha.addActionListener(this);
        nullAlpha.addFocusListener(this);

        btnPanel.add(resetBtn, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
                ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 0, 2, 1), 0, 0));
        btnPanel.add(saveAsBtn, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
                ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 0, 2, 1), 0, 0));

        jPanel1.add(colorscalePanel,  new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 0, 0, 0), 250, 0));
        jPanel1.add(jPanel2,  new GridBagConstraints(0, 3, 1, 1, 1.0, 1.0
            ,GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0));

        if(colorscaleAvailable)
        {
            btnPanel.add(previewBtn,   new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
                ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 0, 2, 1), 0, 0));
            btnPanel.add(cancelBtn, new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0
                ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 0, 2, 1), 0, 0));
        }
        else
            btnPanel.add(cancelBtn, new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
                ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 0, 2, 1), 0, 0));

        paintModePanel.add(paintClr,   new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
                    ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 2, 2, 2), 0, 0));
        paintModePanel.add(indexClr,   new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
                    ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 2, 2, 2), 0, 0));
        paintModePanel.add(activeClr,   new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
                    ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 2, 2, 2), 0, 0));

        modePanel.add(multipleClr, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
                ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 2, 2, 2), 0, 0));
        modePanel.add(keyClr, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
                ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 2, 2, 2), 0, 0));

        jPanel1.add(paintModePanel,  new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
            ,GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0));


        this.getContentPane().add(jPanel1, BorderLayout.WEST);

        btnGroup1.add(multipleClr);
        btnGroup1.add(keyClr);
        btnGroup1.add(interpolateClrs);
        btnGroup2.add(paintClr);
        btnGroup2.add(indexClr);

        jPanel2.add(jLabel1,  new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(6, 5, 5, 0), 8, 6));
        jPanel2.add(jLabel2,  new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(6, 17, 0, 8), 0, 5));
        jPanel2.add(nullRed,  new GridBagConstraints(1, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(6, 13, 5, 0), 19, 1));
        jPanel2.add(nullGreen,  new GridBagConstraints(2, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(6, 8, 5, 0), 19, 1));
        jPanel2.add(jLabel3,  new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(6, 5, 0, 0), 8, 3));
        jPanel2.add(jLabel4,  new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(6, 8, 0, 0), 11, 3));
        jPanel2.add(nullBlue,  new GridBagConstraints(3, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(6, 8, 5, 0), 19, 1));
        jPanel2.add(jLabel5,  new GridBagConstraints(4, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(6, 7, 0, 0), 7, 4));
        jPanel2.add(nullAlpha,  new GridBagConstraints(4, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(6, 7, 5, 0), 19, 1));
        jPanel2.add(nullClrBtn,  new GridBagConstraints(5, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 2, 5), 10, 16));

        jPanel1.add(btnPanel,  new GridBagConstraints(0, 4, 1, 1, 1.0, 1.0
            ,GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 4, 0), 0, 0));
        jPanel1.add(modePanel,  new GridBagConstraints(0, 2, 1, 1, 1.0, 1.0
            ,GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0));
    }

    private void initializeSliderScale()
    {
        maxIndex = nPoints-1;
        minSliderIndex = 0;
        maxSliderIndex = maxIndex;
        if(range != null) {
            minValue = range[0];
            maxValue = range[1];
        }
        minSliderValue = minValue;
        maxSliderValue = maxValue;
        scale = (maxValue - minValue)/maxIndex;
    }

    private void initializeSliderPolygons()
    {
        int totalSliderHeight = sliderHeight+heightPerInterval;
        p = new Point[]
        {
            new Point(-sliderWidth, 0),
            new Point(-sliderWidth, totalSliderHeight),
            new Point( (sliderHeight-sliderWidth), totalSliderHeight),
            new Point(0, heightPerInterval)
        };

        flipX = -1;

        minSliderPolygon = new Polygon();
        minSliderPolygon.addPoint(0, 0);
        for(int n = 0; n < 4; n++)
            minSliderPolygon.addPoint(flipX*p[n].x, p[n].y);

        maxSliderPolygon = new Polygon();
        maxSliderPolygon.addPoint(0, heightPerInterval);
        for(int n = 0; n < 4; n++)
            maxSliderPolygon.addPoint(flipX*p[n].x, -p[n].y+heightPerInterval);

    }

    private void initializeKeyPolygons()
    {
        if(keys == null) return;

        int totalKeyHeight = keyHeight+heightPerInterval;
        p = new Point[]
        {
            new Point( (keyHeight-keyWidth), -totalKeyHeight/2),
            new Point(-keyWidth, -totalKeyHeight/2),
            new Point(-keyWidth, totalKeyHeight/2),
            new Point( (keyHeight-keyWidth), totalKeyHeight/2),
            new Point(0, heightPerInterval)
        };

        flipX = -1;

        keyPolygon = new Polygon[keys.length];
        for(int i=0; i<keys.length; i++) {
            keyPolygon[i] = new Polygon();
            keyPolygon[i].addPoint(0, 0);
            for (int n = 0; n < 5; n++)
                keyPolygon[i].addPoint(flipX * p[n].x, p[n].y);
        }
    }


    private void drawSlider(Graphics g, Polygon sliderPolygon, int sliderIndex, boolean flipVertical)
    {
        String desc = new String();

        int[] px = sliderPolygon.xpoints;
        int[] py = sliderPolygon.ypoints;

        int x = xOrigin;
        int y = yOrigin + barHeight - sliderIndex*heightPerInterval;

        // draw black line from center of bar to slider tip
        g.setColor(Color.black);
        g.drawLine(x, y, x+flipX*barWidth/2, y);

        desc = sliderIndex + " [" + newColors[sliderIndex].getRed() + ", "
            + newColors[sliderIndex].getGreen() + ", "
            + newColors[sliderIndex].getBlue() + ", "
            + newColors[sliderIndex].getAlpha() + "]";
        g.drawString(desc, x+20, y+8);

        g.translate(x, y);

        g.setColor(newColors[sliderIndex]);
        g.fillPolygon(sliderPolygon);

        // highlight above and left, shadow below
        if(flipVertical)
        {
            g.setColor(SystemColor.controlLtHighlight);
            drawHighlight(g, px[1], py[1], px[2], py[2], -flipX, 0);
            drawHighlight(g, px[2], py[2], px[3], py[3], -flipX, -1);
            g.drawLine(px[3], py[3], px[4], py[4]);

            g.setColor(SystemColor.controlDkShadow);
            drawHighlight(g, px[0], py[0], px[1], py[1], 0, 1);
        }
        else
        {
            g.setColor(SystemColor.controlLtHighlight);
            drawHighlight(g, px[0], py[0], px[1], py[1], 0, -1);
            drawHighlight(g, px[1], py[1], px[2], py[2], -flipX, 0);

            g.setColor(SystemColor.controlDkShadow);
            drawHighlight(g, px[2], py[2], px[3], py[3], 0, 1);
            drawHighlight(g, px[3], py[3], px[4], py[4], 0, 1);
        }

        g.translate(-x, -y);
    }

    // draw two lines, one offset by dx,dy from the other
    private void drawHighlight(Graphics g, int x0, int y0, int x1, int y1, int dx, int dy)
    {
        g.drawLine(x0, y0, x1, y1);
        g.drawLine(x0+dx, y0+dy, x1+dx, y1+dy);
    }

    private void drawKeys(Graphics g) {
        int x, y, i;
        String desc = new String();

        int[] px;
        int[] py;

        g.setColor(Color.black);

        x = xOrigin;
        // If there are no keys in the current spectrum create three
        if(keys == null)
        {
            keys = new StsColor[3];
            keys[0] = new StsColor(newColors[0], 0);
            keys[1] = new StsColor(newColors[newColors.length/2], (int)(newColors.length/2));
            keys[2] = new StsColor(newColors[newColors.length - 1], newColors.length - 1);
            initializeKeyPolygons();
        }

        for(i=0; i<keys.length; i++) {
            px = keyPolygon[i].xpoints;
            py = keyPolygon[i].ypoints;
            if(keys[i].idx == -1) {
                int interval = (newColors.length-1)/(keys.length -1);
                keys[i].idx = i * interval;
            }
            if(keys[keys.length-1].idx != newColors.length-1)
                keys[keys.length-1].idx = newColors.length-1;
            y = yOrigin + barHeight - keys[i].idx*heightPerInterval;
            g.drawLine(x, y, x+flipX*barWidth/2+5, y);
            desc = keys[i].idx + " [" + keys[i].getColor().getRed() + ", " + keys[i].getColor().getGreen() + ", " + keys[i].getColor().getBlue() + ", " + keys[i].getColor().getAlpha() + "]";
            if(selectedKeyIndex == i)
                g.drawString("**" + desc, x+20, y+5);
            else
                g.drawString(desc, x+20, y+5);

            g.translate(x, y);

            g.setColor(keys[i].getColor());
            g.fillPolygon(keyPolygon[i]);

            g.setColor(Color.black);
            drawHighlight(g, px[1], py[1], px[2], py[2], -flipX, 0);
            drawHighlight(g, px[2], py[2], px[3], py[3], -flipX, -1);
            g.drawLine(px[3], py[3], px[4], py[4]);

            g.setColor(SystemColor.controlDkShadow);
            drawHighlight(g, px[0], py[0], px[1], py[1], 0, 1);

            g.translate(-x, -y);
        }
    }

    private void addMouseListener()
    {
        MouseListener mouseListener = new MouseAdapter()
        {
            public void mousePressed(MouseEvent e)
            {
                int x, y;
                boolean newIdx = false;

                pickedX = e.getX();
                pickedY = e.getY();
                x = pickedX - xOrigin;
                y = pickedY - (yOrigin + barHeight);
                if((!keyClr.isSelected()) && (minSliderPolygon.contains(x, y + minSliderIndex*heightPerInterval)))
                {
                    sliderSelected = MIN_SLIDER;
                    activeSlider = MIN_SLIDER;
                    pickedIndex = minSliderIndex;
                    startIndex = pickedIndex;
                    activeIndex = pickedIndex;
                    newIdx = true;
                }
                else if((!keyClr.isSelected()) && (maxSliderPolygon.contains(x, y + maxSliderIndex*heightPerInterval)))
                {
                    sliderSelected = MAX_SLIDER;
                    activeSlider = MAX_SLIDER;
                    pickedIndex = maxSliderIndex;
                    activeIndex = pickedIndex;
                    newIdx = true;
                }
                else if(keyClr.isSelected()) {
                    sliderSelected = NONE;
                    activeSlider = NONE;
                    if(keyClr.isSelected()) {
                        for(int i=0; i<keys.length; i++) {
                            if(keyPolygon[i].contains(x, y + keys[i].idx*heightPerInterval)) {
                                selectedKeyIndex = i;
                                newIdx = true;
                            }
                        }
                    if(selectedKeyIndex > -1) pickedIndex = keys[selectedKeyIndex].idx;
                    }
                }
                int mods = e.getModifiers();
                if(newIdx) {
                    if ((mods & InputEvent.BUTTON3_MASK) != 0) {
                        if(sliderSelected == NONE)
                            addKeyPopup(pickedX, pickedY);
                        else
                            addSliderPopup(pickedX, pickedY);
                    }
                } else {
                    pickedIndex = maxIndex - (pickedY - yOrigin)/heightPerInterval;
                    if(pickedIndex > maxIndex)
                        pickedIndex = maxIndex;
                    if(pickedIndex < 0)
                        pickedIndex = 0;
                    activeIndex = pickedIndex;
                    if ((mods & InputEvent.BUTTON3_MASK) != 0) {
                        addColorscalePopup(pickedX, pickedY);
                    } else if(indexClr.isSelected() && !keyClr.isSelected()) {
                        newColors[pickedIndex] = activeColor;
                        saved = false;
                        repaint();
                    }
                }
            }

            public void mouseReleased(MouseEvent e)
            {
                sliderSelected = NONE;
            }
        };
        this.addMouseListener(mouseListener);
    }

    private void addKeyPopup(int x, int y)
    {
        removeListeners();
        JPopupMenu tp = new JPopupMenu("Key Popup");
        colorscalePanel.add(tp);

        tp.add(kc);
        kc.addActionListener(this);
        if((selectedKeyIndex != 0) && (selectedKeyIndex != (keys.length-1))) {
            tp.add(kd);
            kd.addActionListener(this);
        }
        tp.show(colorscalePanel, x, y);
    }

    private void addColorscalePopup(int x, int y) {
        JPopupMenu tp = new JPopupMenu("Colorscale Popup");
        removeListeners();
        colorscalePanel.add(tp);

        if(keyClr.isSelected()) {
            tp.add(ka);
            ka.addActionListener(this);
        }
        tp.add(kr);
        kr.addActionListener(this);
        tp.show(colorscalePanel, x, y);
    }

    private void addSliderPopup(int x, int y) {
        JPopupMenu tp = new JPopupMenu("Slider Popup");
        removeListeners();
        colorscalePanel.add(tp);

        tp.add(kc);
        kc.addActionListener(this);
        tp.show(colorscalePanel, x, y);
    }

    private void removeListeners() {
        kd.removeActionListener(this);
        kc.removeActionListener(this);
        ka.removeActionListener(this);
        kr.removeActionListener(this);
    }

    private void addMouseMotionListener()
    {
        MouseMotionListener motionListener = new MouseMotionAdapter()
        {
            public void mouseDragged(MouseEvent e)
            {
                int i;

                int dy = pickedY - e.getY();
                int index = pickedIndex + dy/heightPerInterval;
                if((index < 255) && (index > -1)) {
                    if (sliderSelected == MIN_SLIDER) {
                        minSliderIndex = index;
                        activeIndex = index;
                        if (minSliderIndex >= maxSliderIndex)
                            minSliderIndex = maxSliderIndex - 1;
                        else if (minSliderIndex < 0)
                            minSliderIndex = 0;
                        computeMinSliderValue();
                    }
                    else if (sliderSelected == MAX_SLIDER) {
                        maxSliderIndex = index;
                        activeIndex = index;
                        if (maxSliderIndex <= minSliderIndex)
                            maxSliderIndex = minSliderIndex + 1;
                        else if (maxSliderIndex > maxIndex)
                            maxSliderIndex = maxIndex;
                        computeMaxSliderValue();
                    }
                    else if (selectedKeyIndex != -1) {

                        // If at an end, set to last index
                        if (index > maxIndex)
                            index = maxIndex;
                        if (index < 0)
                            index = 0;

                            // Must maintain edge keys
                        if (selectedKeyIndex == 0)
                            return;
                        if (selectedKeyIndex == keys.length - 1)
                            return;

                        for (i = selectedKeyIndex + 1; i < keys.length - 1; i++) {
                            if (index > keys[i].idx)
                                keys[i].idx = index;
                        }
                        for (i = 1; i < selectedKeyIndex; i++) {
                            if (index < keys[i].idx)
                                keys[i].idx = index;
                        }
                        keys[selectedKeyIndex].idx = index;

                        // Reinterpolate between keys
                        computeNewColors();
                    }
                    else {
                        if (paintClr.isSelected() && (!keyClr.isSelected())) {
                            int firstIndex = pickedIndex;
                            int lastIndex = index;
                            activeIndex = index;
                            if (firstIndex < 0)
                                firstIndex = 0;
                            if (lastIndex > 254)
                                lastIndex = 254;
                            if (firstIndex > lastIndex) {
                                firstIndex = lastIndex;
                                lastIndex = pickedIndex;
                            }
                            for (i = firstIndex; i < lastIndex; i++) {
                                saved = false;
                                newColors[i] = activeColor;
                            }
                        }
                    }
                    repaint();
                }
            }
        };
        this.addMouseMotionListener(motionListener);
    }

    public void paint(Graphics g)
    {
        super.paint(g);
        g.setFont(graphicFont);

        int nColors = newColors.length - 1;
        barHeight  = nColors * heightPerInterval;

        insets = colorscalePanel.getInsets();
        barRect = colorscalePanel.getBounds();
        barRect.width -= (insets.left + insets.right);
        barRect.height -= (insets.top + insets.bottom);
        barWidth = barRect.width/4;

        heightPerInterval = barRect.height/255;
        xOrigin = barRect.x + barWidth + 10;
        yOrigin = barRect.y + (barRect.height - (heightPerInterval * 255))/2;

        float min = range[0];
        float max = range[1];

        if(nColors > 0) {
            barHeight  = nColors * heightPerInterval;

            // Draw the Colors
            int x = barRect.x + 10;
            int y = yOrigin;

            // Draw Color Scale - Bottom up
            for (int i = nColors - 1; i >= 0; i--) {
                // draw a rectangle for each color
                g.setColor(newColors[nColors - 1 - i]);
                g.fillRect(x, y + (i * heightPerInterval * (255 / nColors)),
                           (int)barWidth, (heightPerInterval * (255 / nColors)));
            }
            if(multipleClr.isSelected()) {
                drawSlider(g, minSliderPolygon, minSliderIndex, false);
                drawSlider(g, maxSliderPolygon, maxSliderIndex, true);
            }

            if(keyClr.isSelected()) {
                drawKeys(g);
            }

            if(multipleClr.isSelected()) {
                g.setColor(Color.black);
                int inty = colorscalePanel.getSize().height/2;
                int intx = colorscalePanel.getSize().width/2;
                String desc = "Index =" + activeIndex;
                g.drawString(desc, intx , inty);
                desc = "Red =" + newColors[activeIndex].getRed();
                g.drawString(desc, intx , inty + 14);
                desc = "Green =" + newColors[activeIndex].getGreen();
                g.drawString(desc, intx , inty + 28);
                desc = "Blue =" + newColors[activeIndex].getBlue();
                g.drawString(desc, intx , inty + 42);
                desc = "Alpha =" + newColors[activeIndex].getAlpha();
                g.drawString(desc, intx , inty + 56);
            }
        }
    }

    public void actionPerformed(ActionEvent e)
    {
        int i;
        Object source = e.getSource();

        // Reset Palette
        if ((source == resetBtn) || (source == kr))
        {
            adjustSpectrum();
            if(origColorscalePanel != null)
            {
                origColorscalePanel.reinitialize();
                model.win3dDisplayAll();
            }
        }
        else if(source == previewBtn) previewSpectrum();
        else if(source == saveAsBtn)  saveSpectrum();
        else if(source == cancelBtn)
        {
            adjustSpectrum();
            if(origColorscalePanel != null)
            {
                origColorscalePanel.reinitialize();
                model.win3dDisplayAll();
            }
            setVisible(false);
        }
        else if(source == activeClr)
        {
            Color color = new Color(255,255,255,255);
            if((pickedIndex > 0) && (pickedIndex < newColors.length))
                color = jColorChoose.showDialog(this,"Color Selection", newColors[pickedIndex]);
            else
                color = jColorChoose.showDialog(this,"Color Selection", color);
            if(color != null)
                activeColor = color;
            activeClr.setBackground(activeColor);
        }
        else if(source == nullClrBtn)
        {
            Color newNullColor = jColorChoose.showDialog(this,"Color Selection", nullColor);
			if(newNullColor == null) return;
			nullColor = newNullColor;
            nullClrBtn.setBackground(nullColor);
            nullRed.setText(new Integer(nullColor.getRed()).toString());
            nullBlue.setText(new Integer(nullColor.getBlue()).toString());
            nullGreen.setText(new Integer(nullColor.getGreen()).toString());
            nullAlpha.setText("255");
        }
        else if((source == nullAlpha) || (source == nullRed) || (source == nullGreen) || (source == nullBlue))
            resetNullValue();
        else if(source == multipleClr)
        {
            // Messqage indicating that changes will be lost
            isOK = true;
            if ((mode == KEY_SELECT) && (saved == false))
                isOK = StsYesNoDialog.questionValue(frame, "You will lose changes made in key mode. Continue?");
            if(!isOK)
            {
                keyClr.setSelected(true);
                return;
            }
            mode = MULTIPLE_SELECT;
            indexClr.setEnabled(true);
            interpolateClrs.setEnabled(true);
            activeClr.setEnabled(true);
            selectedKeyIndex = -1;

        }
        else if(source == keyClr)
        {
            // Messqage indicating that changes will be lost
            isOK = true;
            if(((mode == SINGLE_SELECT) || (mode == MULTIPLE_SELECT)) && (saved == false))
                isOK = StsYesNoDialog.questionValue(frame, "You will lose changes made in Color Mode. Continue?");
            if(!isOK)
            {
                multipleClr.setSelected(true);
                return;
            }
            mode = KEY_SELECT;
            indexClr.setEnabled(false);
            interpolateClrs.setEnabled(false);
            activeClr.setEnabled(false);
            if(keys == null)
                StsMessageFiles.infoMessage("No keys available in this spectrum.");
            else
                computeNewColors();
        }
        else if(source == ka)
        {
            // Add a key
            keys = (StsColor[]) StsMath.arraycopy(keys, keys.length+1);

            //Compute location of new key
            for(i=0; i<keys.length; i++)
            {
                if(keys[i].idx >= pickedIndex)
                {
                    selectedKeyIndex = i - 1;
                    break;
                }
            }
            if(selectedKeyIndex < 0) return;

            for(i = 0; i<=selectedKeyIndex; i++)
                keys[i] = keys[i];

            for(i=keys.length-1; i>=selectedKeyIndex+2; i--)
                keys[i] = keys[i-1];

            keys[selectedKeyIndex + 1] = new StsColor(keys[selectedKeyIndex]);
            keys[selectedKeyIndex + 1].idx = pickedIndex;
            computeNewColors();
            initializeKeyPolygons();

        }
        else if(source == kd)
        {
            // Remove a key
            for(i = selectedKeyIndex; i<keys.length-1; i++)
                keys[i] = keys[i+1];

            keys = (StsColor[]) StsMath.trimArray(keys, keys.length-1);
            computeNewColors();
            initializeKeyPolygons();

        }
        else if(source == kc)
        {
            if(activeSlider == MIN_SLIDER)
            {
                newColors[minSliderIndex] = newColors[minSliderIndex] = jColorChoose.showDialog(this,"Color Selection", newColors[minSliderIndex]);
                interpolateClrs(minSliderIndex,maxSliderIndex);
                saved = false;

            }
            else if (activeSlider == MAX_SLIDER)
            {
                newColors[maxSliderIndex] = jColorChoose.showDialog(this,"Color Selection", newColors[maxSliderIndex]);
                interpolateClrs(minSliderIndex,maxSliderIndex);
                saved = false;

            }
            else if (activeSlider == NONE)
            {
                Color selectedColor = jColorChoose.showDialog(this,"Color Selection", keys[selectedKeyIndex].getColor());
                keys[selectedKeyIndex] = new StsColor(selectedColor);
                computeNewColors();
            }
        }
        Graphics gt = this.getGraphics();
        paint(gt);
        repaint();
    }

    public void focusLost(FocusEvent e)
    {
        Object source = e.getSource();
        if ( (source == nullAlpha) || (source == nullRed) ||
                 (source == nullGreen) || (source == nullBlue))
            resetNullValue();
    }
    public void focusGained(FocusEvent e) { }

    private void resetNullValue()
    {
        Integer tr = new Integer(nullRed.getText());
        int tri = tr.intValue();
        if(tri > 255) tri = 255;
        if(tri < 0) tri = 0;
        Integer tg = new Integer(nullGreen.getText());
        int tgi = tg.intValue();
        if(tgi > 255) tgi = 255;
        if(tgi < 0) tgi = 0;
        Integer tb = new Integer(nullBlue.getText());
        int tbi = tb.intValue();
        if(tbi > 255) tbi = 255;
        if(tbi < 0) tbi = 0;
        Integer ta = new Integer(nullAlpha.getText());
        int tai = ta.intValue();
        if(tai > 255) tai = 255;
        if(tai < 0) tai = 0;

        String val = new Integer(tri).toString();
        nullRed.setText(val);
        val = new Integer(tgi).toString();
        nullGreen.setText(val);
        val = new Integer(tbi).toString();
        nullBlue.setText(val);
        val = new Integer(tai).toString();
        nullAlpha.setText(val);

        nullColor = new Color(tri,tgi,tbi,tai);
        if(spectrum != null) spectrum.setTransparentColor(nullColor);

        nullClrBtn.setBackground(nullColor);
        nullClrBtn.setText("");
        repaint();
    }

    public void previewSpectrum()
    {
        if(colorscaleAvailable)
            colorscale.setNewColors(newColors);
        else
            StsMessageFiles.errorMessage("No object related to this spectrum.");

        origColorscalePanel.repaintModel();
        return;
    }

    public void saveSpectrum()
    {
        String spectrumName = null;

        try
        {
            while(true)
            {
                StsTextAreaDialog dialog = new StsTextAreaDialog(null, "Enter Spectrum Name...", "", 1, 20);
                dialog.setVisible(true);
                spectrumName = dialog.getText();
                if((spectrumName == null) || (spectrumName.length() < 1))
                {
                    StsMessageFiles.infoMessage("Invalid Filename");
                }
                if(model.getObjectWithName(StsSpectrum.class, spectrumName) == null) break;

				if (!StsYesNoDialog.questionValue(model.win3d, "Spectrum name already used, would you like to pick another?")) break;
            }

            StsActionManager windowActionManager = model.mainWindowActionManager;
            StsSpectrum spectrum = null;
            if(!keyClr.isSelected())
            {
				spectrum = StsSpectrum.constructor(spectrumName, newColors);
                if(spectrum == null)
                {
                    StsMessageFiles.errorMessage("Error saving " + spectrumName + " spectrum.");
                    success = false;
                    return;
                }
            }
            else
            {
                spectrum = StsSpectrum.constructor(spectrumName, keys, 255);
            }
            if(colorscale != null)
                colorscale.setSpectrum(spectrum);
            saved = true;
            success = true;
            return;
        }
        catch(Exception e)
        {
            StsException.outputException("StsSpectrumDialog.saveSpectrum() failed.",
                e, StsException.WARNING);
        }
    }

    public void setLabelFormat(String pattern) { labelFormat.applyPattern(pattern); }
    public void setColorscale(StsColorscale colorscale)
    {
        keys = colorscale.getSpectrum().getKeys();
        newColors = (Color []) StsMath.arraycopy(colorscale.getOriginalColors(), colorscale.getNColors());
        this.spectrum = colorscale.getSpectrum();
        this.colorscale = colorscale;
        initialize();
    }
    public void setSpectrum(StsSpectrum spectrum)
    {
        keys = spectrum.getKeys();
        newColors = (Color []) StsMath.arraycopy(spectrum.getColors(), spectrum.getNColors());
        this.spectrum = spectrum;
        this.colorscale = null;
        initialize();
    }

    public int getHeightPerInterval() { return heightPerInterval; }
    public void setHeightPerInterval(int heightPerInterval) { this.heightPerInterval = heightPerInterval; }
    public int getBarWidth() {  return barWidth; }
    public void setBarWidth(int barWidth) { this.barWidth = barWidth; }
    public int getNPoints() { return nPoints; }
    public void setNPoints(int nPoints) { this.nPoints = nPoints; }

    public boolean getSuccess() { return success; }

    public void adjustSpectrum()
    {
        if(spectrum != null)
            newColors = spectrum.getColors();
        else
        {
            newColors[0] = new Color(255,255,255,255);
            newColors[254] = new Color(0,0,0,255);
            interpolateClrs(0, 254);
        }
        repaint();
    }

    private void computeMinSliderValue()
    {
        minSliderValue = minSliderIndex*scale + minValue;
    }

    private void computeMaxSliderValue()
    {
        maxSliderValue = maxSliderIndex*scale + minValue;
    }

    private void interpolateClrs(int first, int last) {
        int nClrs = last - first;
        float[] firstRGBA = new float[4], lastRGBA = new float[4];

        firstRGBA[0] = newColors[first].getRed()/255.0f;
        firstRGBA[1] = newColors[first].getGreen()/255.0f;
        firstRGBA[2] = newColors[first].getBlue()/255.0f;
        firstRGBA[3] = newColors[first].getAlpha()/255.0f;

        lastRGBA[0] = newColors[last].getRed()/255.0f;
        lastRGBA[1] = newColors[last].getGreen()/255.0f;
        lastRGBA[2] = newColors[last].getBlue()/255.0f;
        lastRGBA[3] = newColors[last].getAlpha()/255.0f;
        for(int nn = first; nn <= last; nn++)
        {
            float temp = (float) (nn - first) / (float) nClrs;
            float[] interpolatedColor = StsMath.interpolate(firstRGBA, lastRGBA, temp);
            newColors[nn] = new Color(interpolatedColor[0], interpolatedColor[1],
                                          interpolatedColor[2], interpolatedColor[3]);
        }
    }

    private void computeNewColors() {
        int nColors = newColors.length-1;

        saved = false;
        int nColorIntervals = keys.length - 1;

        int nextColorIndex = 0;
        float[] nextColorRGBA = keys[0].getRGBA();
        float[] scale = new float[3];
        int nn = 0;
        for(int n = 0; n < nColorIntervals; n++)
        {
            int nColorsInInterval = keys[n+1].idx - keys[n].idx;
            int prevColorIndex = nextColorIndex;
            nextColorIndex += nColorsInInterval;
            float[] prevColorRGBA = nextColorRGBA;
            nextColorRGBA = keys[n+1].getRGBA();

            for(nn = prevColorIndex; nn < nextColorIndex; nn++)
            {
                if(nColorsInInterval <= 1) {
                    newColors[nn] = new Color(nextColorRGBA[0], nextColorRGBA[1],
                                              nextColorRGBA[2], nextColorRGBA[3]);
                    nextColorIndex = keys[n+1].idx;
                    nextColorRGBA = keys[n+1].getRGBA();
                    continue;
                } else {
                    float temp = (float) (nn - prevColorIndex) / (float) nColorsInInterval;
                    float[] interpolatedColor = StsMath.interpolate( prevColorRGBA, nextColorRGBA, temp);
                    newColors[nn] = new Color(interpolatedColor[0], interpolatedColor[1],
                                              interpolatedColor[2], interpolatedColor[3]);
                }
            }
        }
        for( ; nn < nColors-1; nn++) {
            newColors[nn] = new Color(nextColorRGBA[0], nextColorRGBA[1],
                                      nextColorRGBA[2], nextColorRGBA[3]);
        }
        if(colorscale != null) colorscale.setNewColors(newColors);
    }

    public void stateChanged(ChangeEvent e)
    {
        Object source = e.getSource();
        return;
    }

    public void mouseExited(MouseEvent e) { }
    public void mouseEntered(MouseEvent e) { }
    public void mousePressed(MouseEvent e) { }
    public void mouseClicked(MouseEvent e) { }
    public void mouseReleased(MouseEvent e) { }
    public void mouseMotion(MouseEvent e) { }

    public void windowClosing(WindowEvent e)
    {
        if(!saved) {
            boolean isOK = StsYesNoDialog.questionValue(frame,
                "Do you want to save this spectrum?");
            if (isOK)
                saveSpectrum();
            else
            {
                // If the preview feature has been used, need to reset to original
                if(colorscaleAvailable)
                    colorscale.setNewColors(colorscale.getOriginalColors());
                dispose();
            }
        }
    }

    public void windowDeactivated(WindowEvent e) { }
    public void windowOpening(WindowEvent e) { }
    public void windowActivated(WindowEvent e) { }
    public void windowDeiconified(WindowEvent e) { }
    public void windowIconified(WindowEvent e) { }
    public void windowClosed(WindowEvent e) { }
    public void windowOpened(WindowEvent e) { }
//
//  Test for dependent constructor
//
    static public void main(String[] args)
    {
        StsColor[] StsColors =
        {
            new StsColor(0, 0, 0, 255, 255),
            new StsColor(255, 0, 255, 255, 255)
        };

        try
        {
            StsModel model = new StsModel();
            StsSpectrum spectrum = StsSpectrum.constructor("Default Palette", StsColors, 255);
            StsColorscale colorscale = new StsColorscale("Test", spectrum, 0.0f, 255.0f);

            StsSpectrumDialog dialog = new StsSpectrumDialog("Test", model, colorscale, true);
            dialog.setVisible(true);
        }
        catch(Exception e)
        {
            e.printStackTrace();
        }
    }
//
//    Test for standalone constructor
//
//    static public void main(String[] args)
//    {
//        StsSpectrumdialog dialog = new StsSpectrumDialog(null, "Test", true);
//        dialog.setVisible(true);
//    }
}
