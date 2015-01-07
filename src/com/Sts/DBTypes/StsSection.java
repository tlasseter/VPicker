//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.Types.Simulation.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.awt.event.*;

public class StsSection extends StsMainObject implements StsSurfaceGridable, StsTreeObjectI, StsTimeDepthDisplayable
{

	// Static members: persistent

	protected int rowDecimation = 4; // not used anymore
	protected int colsPerVertex = 4; // not used anymore

	// Instance members: persistent

	protected StsObjectRefList lines = null; /** List of StsLines defining section: currently two.   */
	protected StsObjectRefList rightLines = null; /** List of lines on right (PLUS) side of section.    */
	protected StsObjectRefList leftLines = null; /** List of lines on left (MINUS) side of section.   */
	protected StsObjectRefList rightBlockSides = null; /** List of blockSides on right (PLUS) side of section.    */
	protected StsObjectRefList leftBlockSides = null; /** List of blockSides on left (MINUS) side of section.   */
	protected StsObjectRefList sectionEdges = null; /** section edges defining geometry of section                 */

//    protected StsSectionEdge sectionEdge = null;        /** edge which defines section at a surface cut (arbitrary for ribbon) */
//    protected StsSectionEdge referenceEdge = null;      /** edge at constant Z which defines section shape          */

	protected StsColor stsColor; /** color of this section (and sectionPatch)    */
	protected float gridUnit = 1.0f; /** Size of a grid cell for defining fault gaps */
	protected StsSectionGeometry geometry;

	// faulted section members
	protected int displaySide = RIGHT; /** current isVisible side relative to + direction */
	protected float[] faultGaps = null; /** faultGaps[0] is RIGHT side and faultGaps[1] is LEFT side */

	// display flags
	protected boolean drawSurface = false;
	protected boolean drawGrid = false;

    /** original domain this object was built in. If and when velocity model is rebuilt, this domain is used as the coordinate source. */
    protected byte zDomainOriginal = StsParameters.TD_NONE;
	// Instance members: transient

	public transient boolean initialized = false; /** Initialization flag */

    transient protected StsList leftZoneSides = null; /**  Zone sides along left side of section                 */
	transient protected StsList rightZoneSides = null; /**  Zone sides along right side of section                 */
    
    /** A faultGrid is constructed for simulation output.  It contains the polygon hierarchy for non-neighbor connections. */
    transient StsFaultGrid faultGrid;

    transient protected int nRows, nCols; /** Size of the section grid (sectionPatch)     */
	transient protected float sectionZInc; /** row Z increment */
	transient protected float sectionZMin = largeFloat; /** min Z value of section */
    transient protected float  sectionZMax = -largeFloat; /** max Z value of section */
    transient StsPatch sectionPatch = null; /** Section grid */
	transient boolean displayOK = true; /** true if drawing this section is OK */
	transient protected StsList gaps = null; /** fault gap info along section */
	transient protected boolean geometryChanged; /** Indicates section geometry has changed: redraw when appropriate */

	transient public boolean drawEdges = true;
	transient public boolean drawLines = true;

	/** Display lists should be used (controlled by View:Display Options) */
	transient boolean useDisplayLists;
	/** Display lists currently being used for surface geometry */
	transient boolean usingDisplayLists = true;

	/** zDomain currently being isVisible. Changing domains requires building new display lists and textures;
	 *  in which case zDomainDisplayed is set to none, display() method deletes display lists, rebuilds displays
	 *  for current project zDomain and sets ZDomainDisplayed to this zDomain.
	 */
	transient protected byte zDomainDisplayed = StsParameters.TD_NONE;

	transient StsBoundingBox boundingBox = null;

	static public StsSectionClass sectionClass = null;
	static protected StsObjectPanel objectPanel = null;
	static protected StsObjectPanel faultObjectPanel = null;

	// convergence parameters used in intersection operations
	static final float maxError = 0.5f;
	static final int maxIter = 20;

	static public final int NONE = StsParameters.NONE;
	static public final int MINUS = StsParameters.MINUS;
	static public final int PLUS = StsParameters.PLUS;

    static public final byte MATCH = StsParameters.MATCH;
	static public final byte MATCH_NOT = StsParameters.MATCH_NOT;
	static public final byte MATCH_REVERSED = StsParameters.MATCH_REVERSED;
	static public final byte MATCH_UNKNOWN = StsParameters.MATCH_UNKNOWN;

//	public static final int VERTICAL_ROW_PLUS = StsParameters.VERTICAL_ROW_PLUS;    /** Vertical ribbon section along row in plus direction.    */
//	public static final int VERTICAL_COL_PLUS = StsParameters.VERTICAL_COL_PLUS;    /** Vertical ribbon section along col in plus direction.    */
//	public static final int VERTICAL_ROW_MINUS = StsParameters.VERTICAL_ROW_MINUS;    /** Vertical ribbon section along row in minus direction.    */
//	public static final int VERTICAL_COL_MINUS = StsParameters.VERTICAL_COL_MINUS;    /** Vertical ribbon section along col in minus direction.    */
//	public static final int VERTICAL_ROW_PLUS_COL_PLUS = StsParameters.VERTICAL_ROW_PLUS_COL_PLUS;    /** Vertical ribbon section along row in plus direction.    */
//	public static final int VERTICAL_ROW_PLUS_COL_MINUS = StsParameters.VERTICAL_ROW_PLUS_COL_MINUS;    /** Vertical ribbon section along col in plus direction.    */
//	public static final int VERTICAL_ROW_MINUS_COL_PLUS = StsParameters.VERTICAL_ROW_MINUS_COL_PLUS;    /** Vertical ribbon section along row in minus direction.    */
//	public static final int VERTICAL_ROW_MINUS_COL_MINUS = StsParameters.VERTICAL_ROW_MINUS_COL_MINUS;    /** Vertical ribbon section along col in minus direction.    */

	/** These are types (type). */
	public static final byte FAULT = StsParameters.FAULT;
	public static final byte BOUNDARY = StsParameters.BOUNDARY;
	public static final byte FRACTURE = StsParameters.FRACTURE;
	public static final byte AUXILIARY = StsParameters.AUXILIARY;

	public static final int BEFORE = StsParameters.BEFORE;
	public static final int AFTER = StsParameters.AFTER;
	public static final int UNKNOWN = StsParameters.UNKNOWN;

	static final int ROW = StsParameters.ROW;
	static final int COL = StsParameters.COL;
	static final int ROWCOL = StsParameters.ROWCOL;

	static final float largeFloat = StsParameters.largeFloat;

	static final int INSIDE = StsParameters.INSIDE;
	static final int OUTSIDE = StsParameters.OUTSIDE;
	static final int NEITHER = 0;

    /** Status flags for type of section; used in constructing simulation grid. */
    static public final byte GEOM_ROW_PLUS = 0;
    static public final byte GEOM_COL_PLUS = 1;
    static public final byte GEOM_ROW_MINUS = 2;
    static public final byte GEOM_COL_MINUS = 3;
    static public final byte GEOM_UNALIGNED = -1;

    private static final float nullValue = StsParameters.nullValue;
	private static final float roundOff = StsParameters.roundOff;
	private static final int nullInteger = StsParameters.nullInteger;

	// display fields
	public static final int RIGHT = StsParameters.RIGHT;
	public static final int LEFT = StsParameters.LEFT;
	static private final String[] sideStrings = {"Left", "Right"};

    static public StsFieldBean[] displayFields = null;
    static public StsFieldBean[] faultDisplayFields = null;
    static public StsFieldBean[] propertyFields = null;
    static public StsFieldBean[] faultPropertyFields = null;

	static protected StsSection currentSection = null; /** current section */

	/** default constructor */
	public StsSection()
	{
        //initializeLists();
	}

	public StsSection(byte type)
	{
		this(type, true);
	}

    public StsSection(boolean persistent)
	{
		super(persistent);
    }

    public StsSection(byte type, boolean persistent)
	{
		super(persistent);
		initializeLists();
		setSectionType(type);
        setZDomain();
	}

	public StsSection(byte type, StsLine firstLine, StsLine lastLine)
	{
		this(type, new StsLine[] { firstLine, lastLine });
	}

    /** use this constructor when you have enough lines to define the columns of the grid (at least 2).
     *  Section may or may not have section edges; if so, they are subsequently added.
     */

    public StsSection(byte type, StsLine[] inputLines)
	{
		this(type, true);
        addLines(inputLines);
        setName(StsParameters.typeToString(type) + "Section-" + getIndex() + " ");
        geometry = new StsSectionGeometry(this, StsSectionGeometry.RIBS, currentModel);
        refreshObjectPanel();
        setZDomain();
	}

    public void addToModel()
	{
		super.addToModel();
		setName(StsParameters.typeToString(type) + "Section-" + getIndex() + " ");
		refreshObjectPanel();
	}

	void initializeLists()
	{
		lines = StsObjectRefList.constructor(2, 2, "lines", this);
		leftLines = StsObjectRefList.constructor(2, 2, "leftLines", this);
		rightLines = StsObjectRefList.constructor(2, 2, "rightLines", this);
	}

    public void addLines(StsLine[] newLines)
    {
        lines.add(newLines);
        ((StsLine)lines.getFirst()).addConnectedSection(this);
        ((StsLine)lines.getLast()).addConnectedSection(this);
    }

    // use this constructor when the section has been completely initialized:
	// lines and edges have been fully defined
	public StsSection(StsSectionEdge edge, byte type, StsLine firstLine, StsLine lastLine, StsSectionGeometry geometry) throws StsException
	{
		this(type, firstLine, lastLine);
		if (edge != null) addSectionEdge(edge);
		this.geometry = geometry;
        setZDomain();
        // completeSection();
		constructSection();
	}

    public StsSection(StsSectionEdge edge, byte type, StsLine firstLine, StsLine lastLine, byte geometryType) throws StsException
	{
		this(type, firstLine, lastLine);
		if (edge != null) addSectionEdge(edge);
        geometry = new StsSectionGeometry(this, geometryType, currentModel);
        setZDomain();
        // completeSection();
		constructSection();
	}

    private void setZDomain()
    {
        if(type != BOUNDARY)
            currentModel.getProject().setModelZDomainToCurrent();
        zDomainOriginal = currentModel.getProject().getZDomain();
    }

    // call this constructor to initially define section: section will be completed later
	static public StsSection constructor(byte type, StsLine firstLine, StsLine lastLine)
	{
		try
		{
			StsSection section = new StsSection(type, firstLine, lastLine);
			return section;
		}
		catch (Exception e)
		{
			StsException.outputException("StsSection.constructorInitializeSection() failed.",
										 e, StsException.WARNING);
			return null;
		}
	}

	// call this constructor when section is fully-defined
	static public StsSection constructor(StsSectionEdge edge, byte type, StsLine firstLine, StsLine lastLine, StsSectionGeometry geometry)
	{
		try
		{
			StsSection section = new StsSection(edge, type, firstLine, lastLine, geometry);
			return section;
		}
		catch (Exception e)
		{
			StsException.outputException("StsSection.constructor() failed.",
										 e, StsException.WARNING);
			return null;
		}
	}

    static public StsSection constructor(StsSectionEdge edge, byte type, StsLine firstLine, StsLine lastLine, byte geometryType)
	{
		try
		{
			StsSection section = new StsSection(edge, type, firstLine, lastLine, geometryType);
			return section;
		}
		catch (Exception e)
		{
			StsException.outputException("StsSection.constructor() failed.",
										 e, StsException.WARNING);
			return null;
		}
	}

    static public StsSection constructor(byte type, StsLine[] lines)
	{
		try
		{
			StsSection section = new StsSection(type, lines);
			return section;
		}
		catch (Exception e)
		{
			StsException.outputException("StsSection.constructor() failed.",
										 e, StsException.WARNING);
			return null;
		}
	}

    public boolean checkConstructSection() throws StsException
	{
		if (sectionPatch != null) return true;
		return completeSection();
		// return constructSection();
	}

	private void computeGeometry()
	{
		if (geometry != null) return;
		geometry = new StsSectionGeometry(this, currentModel);
		geometry.computeDipAngle(this);
	}

    /** section can be drawn if the original domain in which it was constructed is the same as the current zDomain or
     *  a velocityModel exists which can convert from this original zDomain.
     * @return true if it can or can't be drawn for reasons defined above.
     */
    public boolean canDraw()
    {
        if(zDomainOriginal == currentModel.getProject().zDomain) return true;
        return currentModel.getProject().velocityModel != null;
    }

    public boolean isOriginalDepth()
    {
        return zDomainOriginal == StsProject.TD_DEPTH;
    }

    public boolean isZDomainOriginal(byte zDomain)
    {
        return zDomainOriginal == zDomain;
    }

    public boolean constructSection()
	{
        // turn off display of section while computing....
        displayOK = false;
        // set up time or depth range for section
        if(!computeGridRange()) return false;
        // reset edge z values to current domain
//        initializeSectionEdgesTorZ();
        // add points to each edge so they all have points along common columns on the sectionGrid to be constructed
        if(isRibbon())
        {
            constructRibbonSectionEdges();
            constructRibbonSectionGrid();
        }
        else if(isCurved())
        {
            constructCurvedSectionEdges();
            constructCurvedSectionGrid();
        }
        else
        {
            constructRibbedSectionGrid();
            constructRibbedSectionEdges();
        }
        // add interpolated points at each grid crossing so edge follows contour of surface */
        addGridCrossingEdgePoints();
        // set initialized flag and allow for display
        setInitialized(true);
        displayOK = true;
        zDomainDisplayed = currentModel.getProject().getZDomain();
        return true;
	}

    /** called if section edges haven't changed, but zDomain or velocity model has, so rebuild section patch */
/*
    public boolean reconstructSection()
	{
		displayOK = false;
		computeGridRange();
		constructSectionGrid();
//        initializeSectionEdges();
//		if(!geometry.isCurved()) constructRibbonSectionEdges();
		setInitialized(true);
		displayOK = true;
//        zDomainDisplayed = currentModel.project.getZDomain();
		return true;
	}
*/
	// Accessors

	public String getName()
	{
		return super.getName();
	}

	public int getNRows()
	{
		return nRows;
	}

	public int getNCols()
	{
		return nCols;
	}

	public int getRowMin()
	{
		return 0;
	}

	public int getRowMax()
	{
		return nRows - 1;
	}

	public int getColMin()
	{
		return 0;
	}

	public int getColMax()
	{
		return nCols - 1;
	}

	static public StsSection getCurrentSection()
	{
		return currentSection;
	}

	public StsObjectRefList getLines()
	{
		return lines;
	}

	public int getNLines()
	{
		return lines.getSize();
	}

    public StsObjectRefList getSectionEdges()
	{
		if (sectionEdges == null)
		{
			sectionEdges = StsObjectRefList.constructor(2, 2, "sectionEdges", this);
		}
		return sectionEdges;
	}

    public int getNSectionEdges()
    {
        if(sectionEdges == null) return 0;
        return sectionEdges.getSize();
    }

    /*
	 public void setEdge(StsSectionEdge edge)
	 {
	  if (edge == null) return;
	  sectionEdge = edge;
	  if(sectionEdges == null) sectionEdges = StsObjectRefList.constructor(1, 1);
	  sectionEdges.add(edge);
	  edge.setSection(this);
	 }

	 public StsSectionEdge getEdge() { return sectionEdge; }
	 */
	public StsSectionEdge getSectionEdge(StsModelSurface surface)
	{
		if (sectionEdges == null)
		{
			return null;
		}
		int nSectionEdges = sectionEdges.getSize();
		for (int n = 0; n < nSectionEdges; n++)
		{
			StsSectionEdge edge = (StsSectionEdge) sectionEdges.getElement(n);
			if (edge.getSurface() == surface)
			{
				return edge;
			}
		}
		return null;
	}

	public int getDisplaySide()
	{
		return displaySide;
	}

	public StsPatch getPatch()
	{
		return sectionPatch;
	}

	public void setStsColor(StsColor color)
	{
		stsColor = color;
		dbFieldChanged("stsColor", color);
		if (sectionPatch != null)
		{
			sectionPatch.setStsColor(color);
		}
	}

	public StsColor getStsColor()
	{
		return stsColor;
	}

    public boolean isDrawSurface()
	{
		return drawSurface;
	}

	public boolean getDrawSurface()
	{
		return isDrawSurface();
	}

	public void setDrawSurface(boolean drawSurface)
	{
		if (this.drawSurface == drawSurface)
		{
			return;
		}
		this.drawSurface = drawSurface;
		currentModel.win3dDisplayAll();
	}

	public boolean getDrawGrid()
	{
		return drawGrid;
	}

	public void setDrawGrid(boolean drawGrid)
	{
		if (this.drawGrid == drawGrid)
		{
			return;
		}
		this.drawGrid = drawGrid;
		currentModel.win3dDisplayAll();
	}

	public boolean getDrawEdges()
	{
		return drawEdges;
	}

	public void setDrawEdges(boolean drawEdges)
	{
		if (this.drawEdges == drawEdges)
		{
			return;
		}
		this.drawEdges = drawEdges;
		currentModel.win3dDisplayAll();
	}

	/*
	 public void setDrawEdge(Boolean drawEdgesObject)
	 {
	  boolean drawEdges = drawEdgesObject.booleanValue();
	  setDrawEdges(drawEdges);
	 }

	 public void setDrawEdgeIfType(Boolean drawEdgesObject, Integer typeObject)
	 {
	  boolean drawEdges = drawEdgesObject.booleanValue();
	  int type = typeObject.intValue();
	  if(type == this.type) setDrawEdges(drawEdges);
	 }
	 */
	public boolean getDrawLines()
	{
		return drawLines;
	}

	public void setDrawLines(boolean drawLines)
	{
		this.drawLines = drawLines;
		toggleLines(drawLines);
	}

	public StsObjectRefList getLeftBlockSides()
	{
		return leftBlockSides;
	}

	public StsObjectRefList getRightBlockSides()
	{
		return rightBlockSides;
	}

	public StsList getLeftZoneSides()
	{
		return leftZoneSides;
	}

	public StsList getRightZoneSides()
	{
		return rightZoneSides;
	}

    public float[] getXYZorT(int row, int col)
	{
		return sectionPatch.getXYZorT(row, col);
	}

    public float[] getXYZorT(float rowF, float colF)
	{
		return sectionPatch.getXYZorT(rowF, colF);
	}

    public StsPoint getPoint(int row, int col)
	{
		return sectionPatch.getStsPoint(row, col);
	}

	public StsPoint getPoint(float rowF, float colF)
	{
		return sectionPatch.getStsPoint(rowF, colF);
	}

	public void deleteDisplayLists()
	{
		if (sectionPatch == null)
		{
			return;
		}
		GL gl = currentModel.getWin3dGL();
		if (gl == null)
		{
			return;
		}
		sectionPatch.deleteDisplayLists(gl);
	}

	public void setGaps(StsList gaps)
	{
		this.gaps = gaps;
	}

	public StsList getGaps()
	{
		return gaps;
	}

	public StsLine getFirstLine()
	{
		return (StsLine) lines.getFirst();
	}

	public StsLine getLastLine()
	{
		return (StsLine) lines.getLast();
	}

	public StsObjectRefList getRightLines()
	{
		return rightLines;
	}

	public StsObjectRefList getLeftLines()
	{
		return leftLines;
	}

	public boolean isVertical()
	{
        if(geometry == null) return true;
		return geometry.isVertical();
	}

	public float getDipAngle()
	{
		return geometry.getDipAngle();
	}

	public void setInitialized(boolean initialized)
	{
		this.initialized = initialized;
	}

	public boolean getInitialized()
	{
		return initialized;
	}

	public void setGeometry(StsSectionGeometry geometry)
	{
		this.geometry = geometry;
	}

	public StsSectionGeometry getGeometry()
	{
		return geometry;
	}

	public StsLine getOtherEndLine(StsLine line)
	{
		StsLine firstLine = getFirstLine();
		StsLine lastLine = getLastLine();
		if (firstLine == line)
		{
			return lastLine;
		}
		else if (lastLine == line)
		{
			return firstLine;
		}
		else
		{
			return null;
		}
	}

	public boolean hasLines()
	{
		if (leftLines != null && leftLines.getSize() > 0)
		{
			return true;
		}
		if (rightLines != null && rightLines.getSize() > 0)
		{
			return true;
		}
		return false;
	}

	public void setDipAngle(float dipAngle)
	{
		geometry.dipAngleChanged(this, dipAngle);
	}

	// dipAngle has been changed on a number of sections which requires defining
	// lines and connected sections to be updated.  When each section dipAngle
	// was changed, section and the two defining lines were flagged as uninitialized.
	// First for any sections connected to uninitialized lines: flag uninitialized.
	// Second for any sections connected to uninitialized sections: flag uninitialized.

	static public void updateDipAngle()
	{
		int nUninitializedSections;

		boolean initializeIfConnected;
		StsClass sections = currentModel.getCreateStsClass(StsSection.class);

		// flag any sections connected to uninitialized lines as uninitialized
		setSectionsOnUninitializedLines();

		// flag any sections connected to uninitialized Sections as uninitialized
		setSectionsOnUninitializedSections(sections);

		// adjust any lines connected to sections whose dip angle has changed
		adjustLineDipAngles();

		// classInitialize sections not connected to uninitialized sections
		nUninitializedSections = updateDipAngle(sections, initializeIfConnected = false);
		if (nUninitializedSections > 0)
		{
			int noIter = 0;
			while (noIter++ < 3 && nUninitializedSections != 0)
			{
				nUninitializedSections = updateDipAngle(sections, initializeIfConnected = true);
			}
		}

		if (nUninitializedSections > 0)
		{
			StsException.systemError("StsSection.updateDipAngle() failed." +
									 " Unable to classInitialize " + nUninitializedSections + " sections.");
		}
		setNewDipAngle(sections);
	}

	static private void setSectionsOnUninitializedLines()
	{
		StsClass lines = currentModel.getCreateStsClass(StsLine.class);
		int nLines = lines.getSize();
		for (int n = 0; n < nLines; n++)
		{
			StsLine line = (StsLine) lines.getElement(n);
			if (!line.getInitialized())
			{
				StsSection[] connectedSections = line.getAllSections();
				if (connectedSections == null)
				{
					continue;
				}
				for (int s = 0; s < connectedSections.length; s++)
				{
					StsSection section = connectedSections[s];
					section.setInitialized(false);
				}
			}
		}
	}

	static private void setSectionsOnUninitializedSections(StsClass sections)
	{
		// shouldn't have to do this; initialized should already be true...
		/*
		 int nSections = sections.getSize();
		 for(int n = 0; n < nSections; n++)
		 {
		  StsSection section = (StsSection)sections.getElement(n);
		  section.setInitialized(true);
		 }
		 */
		StsClass lines = currentModel.getCreateStsClass(StsLine.class);
		int nLines = lines.getSize();

		for (int n = 0; n < nLines; n++)
		{
			StsLine line = (StsLine) lines.getElement(n);
			StsSection onSection = line.getOnSection();
			if (onSection != null && !onSection.getInitialized())
			{
				line.setInitialized(false);
				StsSection connectedSection = line.getOnlyConnectedSection();
				connectedSection.setInitialized(false);
			}
		}
	}

	static private int updateDipAngle(StsClass sections, boolean initializeIfConnected)
	{
		int nSections = sections.getSize();

		int nUninitializedSections = 0;

		for (int n = 0; n < nSections; n++)
		{
			StsSection section = (StsSection) sections.getElement(n);
			if (!section.updateDipAngle(initializeIfConnected))
			{
				nUninitializedSections++;
			}
		}
		return nUninitializedSections;
	}

	static private void setNewDipAngle(StsClass sections)
	{
		int nSections = sections.getSize();
		for (int n = 0; n < nSections; n++)
		{
			StsSection section = (StsSection) sections.getElement(n);
			section.setNewDipAngle();
		}
	}

	public void setNewDipAngle()
	{
		geometry.setNewDipAngle();
	}

	public float getDipAngleChange()
	{
		return geometry.getDipAngleChange();
	}

	public boolean updateDipAngle(boolean initializeIfConnected)
	{
		StsPoint edgePoint, axis;
		StsLine firstLine, lastLine;
		StsSection onSection;
		boolean firstLineInitialized = true;
		boolean lastLineInitialized = true;

		if (initialized)
		{
			return true;
		}

		try
		{
			float dipAngleChange = getDipAngleChange();

			firstLine = getFirstLine();
			if (!firstLine.getInitialized())
			{
				firstLineInitialized = projectLineToSection(firstLine, initializeIfConnected);

			}
			lastLine = getLastLine();
			if (!lastLine.getInitialized())
			{
				lastLineInitialized = projectLineToSection(lastLine, initializeIfConnected);

			}
			geometry.setIsVertical(false);
			//        setNewDipAngle();

			boolean sectionInitialized = initializeSection();

			boolean linesInitialized = firstLineInitialized && lastLineInitialized;
			initialized = sectionInitialized && linesInitialized;
			return initialized;
		}
		catch (Exception e)
		{
			StsException.outputException("StsSection.updateDipAngle() failed.",
										 e, StsException.WARNING);
			return false;
		}
	}

	static private void adjustLineDipAngles()
	{
		StsClass lines = currentModel.getCreateStsClass(StsLine.class);
		int nLines = lines.getSize();
		for (int w = 0; w < nLines; w++)
		{
			StsLine line = (StsLine) lines.getElement(w);
			if (!line.getInitialized())
			{
				StsSection[] connectedSections = line.getAllSections();
				for (int s = 0; s < connectedSections.length; s++)
				{
					StsSection connectedSection = connectedSections[s];
					float dipAngleChange = connectedSection.getDipAngleChange();
					if (dipAngleChange != 0.0f)
					{
						connectedSection.adjustLineDipAngle(line, dipAngleChange);
					}
				}
			}
		}
	}

	private void adjustLineDipAngle(StsLine line, float dipAngleChange)
	{
		StsPoint edgePoint, axis;
		StsSectionEdge sectionEdge = getPivotSectionEdge();
		if (dipAngleChange != 0.0f)
		{
			if (line == getFirstLine())
			{
				edgePoint = sectionEdge.getFirstPoint();
			}
			else if (line == getLastLine())
			{
				edgePoint = sectionEdge.getLastPoint();
			}
			else
			{
				return;
			}

			axis = getTangentAtPoint(edgePoint, RIGHT);
			line.adjustDipAngle(dipAngleChange, axis.v, edgePoint);
		}
	}

	public StsSectionEdge getSurfaceSectionEdge(StsXYSurfaceGridable surface)
	{
		if (sectionEdges == null)
		{
			return null;
		}
		int nSectionEdges = sectionEdges.getSize();
		for (int n = 0; n < nSectionEdges; n++)
		{
			StsSectionEdge edge = (StsSectionEdge) sectionEdges.getElement(n);
			if (edge.getSurface() == surface)
			{
				return edge;
			}
		}
		return null;
	}

	// return edge about which dip angle is computed; for now make it the first edge picked
	private StsSectionEdge getPivotSectionEdge()
	{
		if (sectionEdges == null)
		{
			return null;
		}
		return (StsSectionEdge) sectionEdges.getFirst();
	}

	public boolean hasSurface(StsXYSurfaceGridable surface)
	{
		if (sectionEdges == null) return false;
		int nSectionEdges = sectionEdges.getSize();
		for (int n = 0; n < nSectionEdges; n++)
		{
			StsSectionEdge edge = (StsSectionEdge) sectionEdges.getElement(n);
			if (edge.getSurface() == surface) return true;
		}
		return false;
	}

	private boolean projectLineToSection(StsLine line, boolean initializeIfConnected)
	{
		StsSection onSection = line.getOnSection();
		if (onSection == null)
		{
			return true;
		}

		boolean onSectionInitialized = onSection.getInitialized();
		if (onSectionInitialized)
		{
			line.projectToSection();
			line.setInitialized(true);
			return true;
		}
		else if (initializeIfConnected)
		{
			line.projectToSection();
			return false;
		}
		else
		{
			return false;
		}
	}

	public float getDZRow()
	{
		return sectionZInc;
	}

	public boolean isBoundary()
	{
		return type == BOUNDARY;
	}

	public boolean isFault()
	{
		return type == FAULT;
	}

	public boolean isAuxiliary()
	{
		return type == AUXILIARY;
	}

	public float getRowCoor(StsGridSectionPoint point)
	{
		return point.getRowF(this);
	}

	public float getRowCoor(StsPoint point)
	{
		return getRowCoor(point.v);
	}

	public float getRowCoor(float[] xyz)
	{
		return getRowF(xyz[2]);
	}

	public float getRowF(float z)
	{
        if(sectionZMin == nullValue) return nullValue;
		float rowF =  (z - sectionZMin) / sectionZInc;
        float newRowF = StsMath.checkRoundOffInteger(rowF);
        if(rowF != newRowF)
            System.out.println("roundoff: rowF " + rowF + " to " + newRowF);
        return newRowF;
    }

	public float getRowValue(float rowF)
	{
		float z = rowF * sectionZInc + sectionZMin;
		return StsMath.minMax(z, sectionZMin, sectionZMax);
	}

	public float getColCoor(StsGridSectionPoint point)
	{
		return point.getColF(this);
	}

	public float getColCoor(StsPoint point)
	{
		return getColF(point);
	}

	public float getColCoor(float[] xyz)
	{
		StsPoint point = new StsPoint(xyz);
		return getColF(point);
	}

	/*
	 public float[] interpolateBilinearXYZ(float[] xyz, boolean computeIfNull)
	 {
	  StsSectionPoint sectionPoint = new StsSectionPoint(new StsPoint(xyz));
	  return sectionPoint.nearestPoint.v;
	 }

	 public float interpolateBilinearZ(float[] xyz, boolean computeIfNull)
	 {
	  StsSectionPoint sectionPoint = new StsSectionPoint(new StsPoint(xyz));
	  return sectionPoint.nearestPoint.v[2];
	 }
	 */
	/*
	 public float getColCoor(StsPoint point)
	 {
	  float f;

	  float[] xyzf = point.v;
	  if(xyzf == null) return nullValue;

	  if(xyzf.length < 4)
	   f = getSectionF(point);
	  else
	   f = xyzf[3];

	  return getColF(f);
	 }
	 */
	/*
	 public float getColF(float f)
	 {
	  if(f < 0.0f || f > 1.0f) return nullValue;
	  return f*(nCols - 1);
	 }
	 */
// These methods don't fit interface: redefine interface?

	public float getColValue(float colF)
	{
		return colF;
	}

	public float getXInc()
	{
		return nullValue;
	}

	public float getYInc()
	{
		return nullValue;
	}

	public float getXMax()
	{
		return nullValue;
	}

	public float getXMin()
	{
		return nullValue;
	}

	public float getYMax()
	{
		return nullValue;
	}

	public float getYMin()
	{
		return nullValue;
	}

	public float getZMax()
	{
		return sectionZMax;
	}

	public float getZMin()
	{
		return sectionZMin;
	}

	public float getZAtRow(int row)
	{
		return sectionZMin + sectionZInc * row;
	}

	public StsObjectRefList getLinesOnSectionSide(int side)
	{
		if (side == LEFT)
		{
			return leftLines;
		}
		else if (side == RIGHT)
		{
			return rightLines;
		}
		else
		{
			return null;
		}
	}

	private void toggleLines(boolean on)
	{
		int nLines = (lines == null) ? 0 : lines.getSize();
		for (int i = 0; i < nLines; i++)
		{
			StsLine line = (StsLine) lines.getElement(i);
			line.setIsVisible(on);
		}
	}

	public boolean isRibbon()
	{
		return geometry.isRibbon();
	}

	public boolean isCurved()
	{
		return geometry.isCurved();
	}

    public boolean isRibbed()
    {
        return geometry.isRibbed();
    }

	public boolean isPlanar()
	{
		return geometry.isPlanar();
	}

	public boolean isRowColAligned()
	{
		return geometry.isRowColAligned();
	}

    public boolean edgesAreStraight()
    {
        int nSectionEdges = sectionEdges.getSize();
        for(int n = 0; n < nSectionEdges; n++)
        {
            StsSectionEdge edge = (StsSectionEdge)sectionEdges.getElement(n);
            int nEdgePoints = edge.edgePoints.getSize();
            if(nEdgePoints > 2) return false;
        }
        return true;
    }

    // Section can be vertical or non-vertical.  If vertical return whether it is
	// aligned along ROW or COL.  If vertical, but not so aligned, return ROWCOL.
	// If nonVertical, return NONE.

	public int getRowColAlignment()
	{
		if (geometry == null)
		{
			return NONE;
		}
		else
		{
			return geometry.getRowOrCol();
		}
	}

	public StsLine getSectionLine(int end)
	{
		if (end == MINUS)
		{
			return getFirstLine();
		}
		else if (end == PLUS)
		{
			return getLastLine();
		}
		else
		{
			StsException.outputException(new StsException(StsException.WARNING,
				"StsSection.getSectionLine(int)", " called with wrong argument: " + end +
				"; must be MINUS(-1) or PLUS(+1)."));
			return null;
		}
	}

	public float[] getPlanarNormal()
	{
		return sectionPatch.getPlanarNormal();
	}

	public float[] getNormal(int row, int col)
	{
		return sectionPatch.getNormal(row, col);
	}

	public float[] getNormal(float rowF, float colF)
	{
		return sectionPatch.getNormal(rowF, colF);
	}

	/*
	 public void construct()
	 {
	  try
	  {
	   if (isRibbon()) constructRibbonSection();
	   else constructCurvedSection();
//			constructIntersectedEdges();
	  }
	  catch(Exception e)
	  {
	   StsException.outputException("StsSection.construct() failed.",
		e, StsException.WARNING);
	  }
	 }
	 */

	public void setRightGap(float gap)
	{
		if (faultGaps == null)
		{
			initFaultGaps();
		}
		if (faultGaps[0] == gap)
		{
			return;
		}
		faultGaps[0] = gap;
		dbFieldChanged("faultGaps", faultGaps);
	}

	public void setLeftGap(float gap)
	{
		if (faultGaps == null)
		{
			initFaultGaps();
		}
		if (faultGaps[1] == gap)
		{
			return;
		}
		faultGaps[1] = gap;
		dbFieldChanged("faultGaps", faultGaps);
	}

	public float getRightGap()
	{
		if (faultGaps == null)
		{
			initFaultGaps();
		}
		return faultGaps[0];
	}

	public float getLeftGap()
	{
		if (faultGaps == null)
		{
			initFaultGaps();
		}
		return faultGaps[1];
	}

	public float[] getFaultGaps()
	{
		if (faultGaps == null)
		{
			initFaultGaps();
		}
		return faultGaps;
	}

	private void initFaultGaps()
	{
		float gridGaps;

		faultGaps = new float[2];

		if (type == FAULT)
		{
			gridGaps = 0.5f;
		}
		else
		{
			gridGaps = 0.0f;

		}
		faultGaps[0] = gridGaps;
		faultGaps[1] = gridGaps;
	}

	public float getFaultGap(int side)
	{
		if (faultGaps == null)
		{
			initFaultGaps();
		}
		if (side == RIGHT)
		{
			return faultGaps[0];
		}
		else
		{
			return faultGaps[1];
		}
	}

	/*
	 public float[] getFaultGapDistances(StsSurface surface)
	 {
	  gapDistances = new float[2];
	  float xInc, yInc;

	  if( surface != null )
	  {
	   StsGrid grid = surface.getGrid();
	   xInc = grid.getXInc();
	   yInc = grid.getYInc();

	   if(xInc == yInc)
		gridUnit = xInc;
	   else
		gridUnit = (float)Math.sqrt(xInc*yInc);
	  }
	  gapDistances[0] = faultGaps[0]*gridUnit;
	  gapDistances[1] = faultGaps[1]*gridUnit;
	 }
	 */


	public boolean initialize(StsModel model)
	{
//		return true;
		return initialize();
	}

	public boolean initialize()
	{
		try
		{
			if (initialized) return true;

            if(sectionClass == null) sectionClass = (StsSectionClass)currentModel.getStsClass(StsSection.class);
            lines.forEach("initialize");

            StsObject[] lineObjects = lines.getElements();
            boolean linesOk = true;
            for(int n = 0; n < 2; n++)
            {
                StsLine line = (StsLine)lineObjects[n];
                if(!line.initialize()) linesOk = false;
            }
            if(!linesOk) return false;
            initialized = initializeSection();
            return initialized;
        }
		catch (Exception e)
		{
			StsException.outputException("Failed to classInitialize section " + getName(), e, StsException.WARNING);
			delete();
			initialized = false;
            return initialized;
        }
	}

	// recompute sectionGrid; if endSections are initialized, then this section is initialized
	// return true if successfully completed (though maybe not initialized)
	public boolean reinitialize()
	{
		try
		{
			displayOK = false;
			if (!initializeSectionGrid())
			{
				return false;
			}
			displayOK = true;
			// set initialized true if end sections are already initialized
			return isEndSectionsInitialized();
		}
		catch (Exception e)
		{
			StsException.outputException("Failed to classInitialize section", e, StsException.WARNING);
			initialized = false;
			return false;
		}
	}

	public boolean isEndSectionsInitialized()
	{
		StsSection endSection;
		endSection = getFirstLine().getOnSection();
		if (endSection != null && !endSection.initialized)
		{
			return false;
		}
		endSection = getLastLine().getOnSection();
		if (endSection != null && !endSection.initialized)
		{
			return false;
		}
		return true;
	}

	public void checkDelete()
	{
		if (lines == null || lines.getFirst() == null || lines.getLast() == null)
        {
            StsMessageFiles.errorMessage("section: " + getName() + " deleted because terminating fault lines are missing.");
            delete();
        }
    }

	public boolean checkSection(StsWin3d win3d)
	{
		if (lines.getSize() < 2)
		{
//            StsTransaction.currentAbort();
			new StsMessage(win3d, StsMessage.WARNING, "Less than 2 lines picked\n" +
						   "for the section:\n" + "construction aborted.");
			return false;
		}
		if (duplicateSection())
		{
//            StsTransaction.currentAbort();
			new StsMessage(currentModel.win3d, StsMessage.WARNING, "A section already has\n" +
						   "these two line/faults:\n" + "construction aborted.");
			return false;
		}

		return true;
	}

	public boolean hasAllModelSurfaceEdges()
	{
		StsModelSurface[] modelSurfaces = (StsModelSurface[]) currentModel.getCastObjectList(StsModelSurface.class);
		if (modelSurfaces == null)
		{
			return false;
		}
		if (sectionEdges == null)
		{
			return false;
		}
		return sectionEdges.getSize() == modelSurfaces.length;
	}

	public boolean completeSection()
	{
		setDisplay();
		initFaultGaps();
		computeGeometry();
		setLines();
		geometry.computeDipAngle(this);
		setCurrentSection();
		if (currentModel != null)
		{
			currentModel.getSpectrumClass().incrementSpectrumColor("Basic");
		}
        return constructSection();
    }

	public void addLineChangesToDB()
	{
		for (int n = 0; n < lines.getSize(); n++)
		{
			StsLine line = (StsLine) lines.getElement(n);
			currentModel.instanceChange(line, line.getName() + " changed");
		}
	}

	public float computeDipAngle(float[] dipDirectionVector)
	{
		StsSectionEdge sectionEdge = getPivotSectionEdge();
        if(sectionEdge == null) return 0.0f;
        StsPoint[] edgePoints = sectionEdge.getPoints();
		if (edgePoints == null) return 0.0f;

		int nEdgePoints = edgePoints.length;

		float dipAngle1 = getFirstLine().getDipAngle(dipDirectionVector, edgePoints[0], MINUS);
		float dipAngle2 = getLastLine().getDipAngle(dipDirectionVector, edgePoints[nEdgePoints - 1], PLUS);
		return 0.5f * (dipAngle1 + dipAngle2);
	}

	public StsSection duplicateSection(StsSection[] sections)
	{
		if (sections == null)
		{
			return null;
		}
		int nSections = sections.length;
		for (int n = 0; n < nSections; n++)
		{
			StsSection duplicateSection = sections[n];
			if (duplicateLines(duplicateSection))
			{
				return duplicateSection;
			}
		}
		return null;
	}

	private boolean duplicateSection()
	{
		if (currentModel == null)
		{
			return false;
		}
		StsClass sections = currentModel.getCreateStsClass(StsSection.class);

		StsSection section = (StsSection) sections.getFirst();

		while (section != null)
		{
			if (section != this && duplicateLines(section))
			{
				return true;
			}
			section = (StsSection) sections.getNext();
		}

		return false;
	}

	private boolean duplicateLines(StsSection section)
	{
		StsLine[] otherLines = section.getLineList();
		StsLine[] thisLines = this.getLineList();

		if (otherLines[0] == thisLines[0] && otherLines[1] == thisLines[1] ||
			otherLines[0] == thisLines[1] && otherLines[1] == thisLines[0])
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	private void setSectionType(byte type)
	{
		this.type = type;
		setName(StsParameters.typeToString(type) + "Section-" + getIndex() + " ");
		if(type == FAULT) setStsColor(StsSectionClass.defaultFaultColor);
	}

	public void setLines()
	{
		addLineToSectionSide(MINUS);
		addLineToSectionSide(PLUS);
	}

	public void setLine(StsLine line)
	{
		int end = getEnd(line);
		if (end == NONE)
		{
			return;
		}
		addLineToSectionSide(end);
	}

	public void setLineOnSide(StsLine line)
	{
		StsSection connectedSection = line.getOnlyConnectedSection();
		int end = connectedSection.getEnd(line);
		if (end == NONE)
		{
			return;
		}
		connectedSection.addLineToSectionSide(end);
	}

	public void addLineToSectionSide(int end)
	{
		StsLine line;

		StsSectionEdge sectionEdge = getPivotSectionEdge();
		line = getSectionLine(end);
		StsSection onSection = line.getOnSection();
		if (onSection == null)
		{
			return;
		}

		StsPoint sectionTangent = sectionEdge.computeEndPointTangent(end);
		StsPoint commonPoint = getEdgeEndPoint(sectionEdge, end);
		StsPoint oldSectionTangent = getTangentAtPoint(commonPoint);
		StsPoint onSectionTangent = onSection.getTangentAtPoint(commonPoint);

		float side = StsPoint.crossProduct2D(onSectionTangent, sectionTangent);
		if (end == MINUS && side < 0.0f || end == PLUS && side > 0.0f)
		{
			addLineToSectionSide(line, end, onSection, RIGHT);
		}
		else
		{
			addLineToSectionSide(line, end, onSection, LEFT);
		}
	}

	private boolean addLineToSectionSide(StsLine line, int end, StsSection onSection, int side)
	{
		StsObjectRefList sideLines;
		StsLine prevLine = null, nextLine;

		line.setSectionSide(side);

		sideLines = onSection.getLinesOnSectionSide(side);
		int nLines = sideLines.getSize();

		if (nLines == 0)
		{
			sideLines.add(line);
			return true;
		}

		// eliminate any bad sideLines:  this is a hacque
		for (int n = 0; n < nLines; n++)
		{
			StsLine sectionLine = (StsLine) sideLines.getElement(n);
			if (sectionLine == null || sectionLine.getIndex() < 0)
			{
				sideLines.delete(n);
			}
		}

		/** Check if line has already been added */
		for (int n = 0; n < nLines; n++)
		{
			if (line == sideLines.getElement(n))
			{
				return false;
			}
		}

		nextLine = null;
		for (int n = 0; n < nLines; n++)
		{
			prevLine = nextLine;
			nextLine = (StsLine) sideLines.getElement(n);
			int position = compareLinePositions(line, nextLine, onSection, side);

			if (position == UNKNOWN)
			{
				StsException.outputException(new StsException(StsException.WARNING,
					"StsSection.addLineToSectionSide() failed.", "Position of side edge " +
					line.getIndex() + " relative to edge " + nextLine.getIndex() +
					" can't be determined."));
				return false;
			}
			else if (position == BEFORE)
			{
				if (!checkAdjustLineVertexPositions(line, prevLine, nextLine, onSection, side))
				{
					return false;
				}
				sideLines.insertBefore(nextLine, line);
				return true;
			}
		}

		prevLine = nextLine;
		if (!checkAdjustLineVertexPositions(line, prevLine, null, onSection, side))
		{
			return false;
		}
		sideLines.add(line);
		return true;
	}

	/** return position of lineA relative to lineB on section
	 *  by comparing sectionCol values. Check vertices on surfaces only.
	 *  On RIGHT side, col is from 0 to nCols-1 and on LEFT, col is from nCols-1 to 0.
	 */
	private int compareLinePositions(StsLine lineA, StsLine lineB, StsSection section, int side)
	{
		int position, newPosition;
		StsObjectRefList lineVertices;
		int nLineVertices;
		boolean extrapolate = true;

		position = UNKNOWN;

		lineVertices = lineA.getSectionEdgeVertices();
		if (lineVertices == null)
		{
			lineVertices = lineA.getLineVertices();
		}
		nLineVertices = lineVertices.getSize();

		int nBefore = 0;
		int nAfter = 0;
		for (int n = 0; n < nLineVertices; n++)
		{
			StsSurfaceVertex vertex = (StsSurfaceVertex) lineVertices.getElement(n);

			float colA = vertex.getColF(section);
			float aZ = vertex.getPoint().getZorT();
			StsPoint bPnt = lineB.getXYZPointAtZorT(aZ, extrapolate);
			float colB = section.getColF(bPnt);

			if (colA == colB)
			{
				continue;
			}
			if (side == RIGHT)
			{
				newPosition = (colA < colB) ? BEFORE : AFTER;
			}
			else
			{
				newPosition = (colA > colB) ? BEFORE : AFTER;

			}
			if (newPosition == BEFORE)
			{
				nBefore++;
			}
			if (newPosition == AFTER)
			{
				nAfter++;
			}
		}

		lineVertices = lineB.getLineVertices();
		nLineVertices = lineVertices.getSize();
		for (int n = 0; n < nLineVertices; n++)
		{
			StsSurfaceVertex vertex = (StsSurfaceVertex) lineVertices.getElement(n);
			/*
			 StsSurface surface = vertex.getSurface();
			 if (surface == null)
			 {
			  continue;
			 }
			 */
			float colB = vertex.getColF(section);
			float bZ = vertex.getPoint().getZorT();
			StsPoint aPnt = lineA.getXYZPointAtZorT(bZ, extrapolate);
			float colA = section.getColF(aPnt);

			if (side == RIGHT)
			{
				newPosition = (colA < colB) ? BEFORE : AFTER;
			}
			else
			{
				newPosition = (colA > colB) ? BEFORE : AFTER;

			}
			if (newPosition == BEFORE)
			{
				nBefore++;
			}
			if (newPosition == AFTER)
			{
				nAfter++;
			}
		}
		if (nBefore > nAfter)
		{
			return BEFORE;
		}
		else if (nBefore < nAfter)
		{
			return AFTER;
		}
		else
		{
			return UNKNOWN;
		}
	}

	// Line has been checked for surface vertices in right position relative to
	// prevLine and nextLine on this side.  Now adjust the other vertices so they
	// don't cross over neighboring lines on side. If they do, adjust so they are
	// between limits.  If prevLine exists, one limit is col of prevLine at same z
	// Return true if operation performed without error; false otherwise.
	private boolean checkAdjustLineVertexPositions(StsLine line,
		StsLine prevLine, StsLine nextLine,
		StsSection section, int side)
	{
		int position, newPosition;
		StsObjectRefList lineVertices;
		int nLineVertices;
		boolean extrapolate = true;
		float colFLimit = 0.0f, colF, prevColF, nextColF;
		int nTopVertex, nBotVertex;
		StsSurfaceVertex vertex;
		StsXYSurfaceGridable surface;
		StsPoint point;
		boolean adjusted = false;

		try
		{
			position = UNKNOWN;

			lineVertices = line.getLineVertices();
			nLineVertices = lineVertices.getSize();

			// adjust points above top vertex on a surface
			nTopVertex = 0;
			for (int n = 0; n < nLineVertices; n++)
			{
				vertex = (StsSurfaceVertex) lineVertices.getElement(n);
				surface = vertex.getSurface();
				if (surface == null)
				{
					continue;
				}
				colFLimit = vertex.getColF(section);
				nTopVertex = n;
				break;
			}
			for (int n = 0; n < nTopVertex; n++)
			{
				vertex = (StsSurfaceVertex) lineVertices.getElement(n);
				if (adjustLineVertex(vertex, prevLine, nextLine, section, side, colFLimit))
				{
					adjusted = true;
				}
			}

			// adjust points below bot vertex on a surface
			nBotVertex = nLineVertices - 1;
			for (int n = nLineVertices - 1; n >= 0; n--)
			{
				vertex = (StsSurfaceVertex) lineVertices.getElement(n);
				surface = vertex.getSurface();
				if (surface == null)
				{
					continue;
				}
				colFLimit = vertex.getColF(section);
				nBotVertex = n;
				break;
			}
			for (int n = nBotVertex + 1; n < nLineVertices; n++)
			{
				vertex = (StsSurfaceVertex) lineVertices.getElement(n);
				if (adjustLineVertex(vertex, prevLine, nextLine, section, side, colFLimit))
				{
					adjusted = true;
				}
			}
			// there may be vertices in between surfaces:  add checks for those

			if (adjusted)
			{
				line.computePointsProjectToSection();
			}
			return true;
		}
		catch (Exception e)
		{
			StsException.outputException("StsSection.checkAdjustLineVertexPositions() failed.",
										 e, StsException.WARNING);
			return false;
		}
	}

	private boolean adjustLineVertex(StsSurfaceVertex vertex, StsLine prevLine,
									 StsLine nextLine, StsSection section, int side, float colFLimit)
	{
		StsPoint point;
		float colF, prevColF, nextColF, newColF;
		float z;
		boolean extrapolate = true;
		boolean adjusted = false;

		try
		{
			colF = vertex.getColF(section);
			z = vertex.getPoint().getZorT();
			prevColF = getLineColForSide(prevLine, z, section, side, MINUS);
			nextColF = getLineColForSide(nextLine, z, section, side, PLUS);

			if (StsMath.betweenInclusive(colF, prevColF, nextColF))
			{
				return adjusted;
			}

			if (prevColF > nextColF)
			{
				float tempColF = prevColF;
				prevColF = nextColF;
				nextColF = tempColF;
			}
			if (colF <= prevColF)
			{
				float maxColF = Math.min(nextColF, colFLimit);
				newColF = 0.5f * (prevColF + maxColF);
				adjustVertexColF(vertex, newColF, section);
				adjusted = true;
			}
			else if (colF >= nextColF)
			{
				float minColF = Math.max(prevColF, colFLimit);
				newColF = 0.5f * (nextColF + minColF);
				adjustVertexColF(vertex, newColF, section);
				adjusted = true;
			}
			return adjusted;
		}
		catch (Exception e)
		{
			StsException.outputException("StsSection.adjustLineVertex() failed.",
										 e, StsException.WARNING);
			return adjusted;
		}
	}

	private float getLineColForSide(StsLine line, float z, StsSection section, int side, int direction)
	{
		if (line != null)
		{
			StsPoint point = line.getXYZPointAtZorT(z, true);
			return section.getColF(point);
		}
		else if (side == RIGHT)
		{
			if (direction == MINUS)
			{
				return 0.0f;
			}
			else
			{
				return section.getColMax();
			}
		}
		else
		{
			if (direction == PLUS)
			{
				return 0.0f;
			}
			else
			{
				return section.getColMax();
			}
		}
	}

	private void adjustVertexColF(StsSurfaceVertex vertex, float newColF, StsSection section)
	{
		float rowF = vertex.getRowF(section);
		StsPoint point = section.getPoint(rowF, newColF);
		vertex.getSurfacePoint().setPoint(point);
	}

	/** return position of lineA relative to lineB on section
	 *  by comparing sectionCol values. Check all vertices for consistency.
	 *  On RIGHT side, f is from 0 to 1 and on LEFT, f is from 1 to 0.
	 *  Check only between zRange shared by both line/faults.
	 */

	/** Set the initial display side to RIGHT; fault sections can be toggled;
	 *  other sides have only one side: arbitrarily the RIGHT.
	 */
	private void setDisplay()
	{
		displaySide = RIGHT;
		if (type == FAULT)
		{
			drawSurface = true;
		}
	}

	public void setFaultGaps(float rightGap, float leftGap)
	{
		if (faultGaps == null)
		{
			initFaultGaps();
		}
		faultGaps[0] = rightGap;
		faultGaps[1] = leftGap;
	}

	public void setDisplaySide(int displaySide)
	{
		this.displaySide = displaySide;
	}

	public void setDisplaySideName(String sideName)
	{
		if (sideName.equals("Left"))
		{
			displaySide = LEFT;
		}
		if (sideName.equals("Right"))
		{
			displaySide = RIGHT;
		}
		dbFieldChanged("displaySide", displaySide);
	}

	public String getDisplaySideName()
	{
		if (displaySide == LEFT)
		{
			return "Left";
		}
		else
		{
			return "Right";
		}
	}

	public void addSectionsToLines()
	{
		lines.forEach("addConnectedSection", this);
	}

	public void toggleDrawSurface(ItemEvent e)
	{
		drawSurface = (e.getStateChange() == ItemEvent.SELECTED);
		try
		{
			currentModel.win3d.win3dDisplay();
		}
		catch (NullPointerException ex)
		{}
	}

	public void showGapsDialog(StsWin3d win3d)
	{
		// the dialog displays the gaps in grid units
		StsSectionGapDialog d = new StsSectionGapDialog();
		d.setRightGap(getFaultGap(RIGHT));
		d.setLeftGap(getFaultGap(LEFT));
		d.setGapColor(stsColor.getColor());

		d.setLocationRelativeTo(win3d);
		d.setVisible(true);
		if (d.okWasPressed())
		{
			System.out.println("Setting Gap values> Right: " +
							   d.getRightGap() * gridUnit + " | left: " + d.getLeftGap() * gridUnit);
			setFaultGaps(d.getRightGap() * gridUnit, d.getLeftGap() * gridUnit);
		}
		d.dispose();
	}

	public boolean initializeSection()
	{
        /*
        if (sectionEdges != null)
		{
			sectionEdges.forEach("constructPoints", null);
		}
        */
		return constructSection();
	}

	private boolean computeGridRange()
	{
		StsPoint[][] sectionPoints;
		StsLine[] sectionLines;
		int n;

		try
		{
			sectionZMin = largeFloat;
			sectionZMax = -largeFloat;

			sectionLines = getLineList();
			if (sectionLines == null)
			{
				return false;
			}

			for (n = 0; n < 2; n++)
			{
                StsLine sectionLine = sectionLines[n];
                if(sectionLine == null) return false;
                if(!sectionLine.initialized) sectionLine.initialize();
                float topZ = sectionLine.getTopPoint().getZorT();
				float botZ = sectionLine.getBotPoint().getZorT();
//                StsPoint[] points = sectionLines[n].getPoints();
                // System.out.println("sectionZMin "+sectionZMin + " topZ " + topZ);
				if(topZ < sectionZMin) sectionZMin = topZ;
                // System.out.println("sectionZMin "+ sectionZMin);
				sectionZMax = Math.max(sectionZMax, botZ);
			}

			if (sectionEdges != null)
			{
				int nSectionEdges = sectionEdges.getSize();
				for (n = 0; n < nSectionEdges; n++)
				{
					StsSectionEdge edge = (StsSectionEdge) sectionEdges.getElement(n);
                    float[] zRange = edge.getRangeZorT();
                    if(zRange[0] < sectionZMin) sectionZMin = zRange[0];
                    if(zRange[1] > sectionZMax) sectionZMax = zRange[1];
                    StsXYSurfaceGridable surface = edge.getSurface();
                    if(surface != null)
                    {
	                    float zMin = surface.getZMin();
	                    if(zMin != nullValue)
	                        sectionZMin = Math.min(sectionZMin, zMin);
						float zMax = surface.getZMax();
	                    if(zMax != nullValue)
	                        sectionZMax = Math.max(sectionZMax, zMax);
                    }
				}
            }
            
            StsProject project = currentModel.getProject();

			int sectionIndexMin = project.getIndexAbove(sectionZMin, isDepth);
			int sectionIndexMax = project.getIndexBelow(sectionZMax, isDepth);

			// Reset sectionZMin and sectionZMax to correspond to rowDecimated limits
			sectionZMin = project.getZAtIndex(sectionIndexMin, isDepth);
			sectionZMax = project.getZAtIndex(sectionIndexMax, isDepth);
//            dbFieldChanged("sectionZMin", sectionZMin);
//            dbFieldChanged("sectionZMax", sectionZMax);
            if (isVertical())
			{
				nRows = 2;
				sectionZInc = sectionZMax - sectionZMin;
			}
			else
			{
				nRows = (sectionIndexMax - sectionIndexMin + 1);
				sectionZInc = project.getZorTInc(isDepth);
			}
//            dbFieldChanged("nRows", nRows);
//            dbFieldChanged("dZRow", dZRow);
            return true;
        }
		catch (Exception e)
		{
			StsException.outputWarningException(this, "computeGridRange", e);
            return false;
        }
	}

	public StsLine[] getLineList()
	{
		StsLine[] sectionLines = new StsLine[2];
		sectionLines[0] = (StsLine) lines.getFirst();
//        if(sectionLines[0] == null) return null;
		sectionLines[1] = (StsLine) lines.getLast();
//        if(sectionLines[1] == null) return null;
		return sectionLines;
	}

	public StsPoint[] getSectionSidePointsAtZ(float z, boolean extrapolate)
	{
		try
		{
			StsLine[] sectionLines = getLineList();
			if (sectionLines == null)
			{
				return null;
			}

			StsPoint[] sectionSidePoints = new StsPoint[2];

			for (int n = 0; n < 2; n++)
			{
				sectionSidePoints[n] = sectionLines[n].getXYZPointAtZorT(z, extrapolate);

			}
			return sectionSidePoints;
		}
		catch (Exception e)
		{
			StsException.outputException("StsSeciton.getSectionSidePointsAtZ() failed.",
										 e, StsException.WARNING);
			return null;
		}
	}

	// normals are unconditionally built for sections
	public void checkConstructGridNormals()
	{
		sectionPatch.checkMakeNormals();
	}

	public void constructRibbonSectionEdges()
	{
		try
        {
            StsModelSurface[] surfaces = (StsModelSurface[]) currentModel.getCastObjectList(StsModelSurface.class);
            int nSurfaces = surfaces.length;
            if (nSurfaces == 0)
            {
                return;
            }

            if (sectionEdges == null)
            {
                sectionEdges = StsObjectRefList.constructor(nSurfaces, 1, "sectionEdges", this);
            }
            for (int n = 0; n < nSurfaces; n++)
            {
                if (hasSurface(surfaces[n])) continue;
                StsSectionEdge edge = new StsSectionEdge(type, this, surfaces[n], 1);
                sectionEdges.add(edge);
            }
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "constructRibbonSectionEdges", e);
        }
    }

	/** Build a curved section from a single picked edge.
	 * If this is already a ribbon section, delete the
	 * single ribbon section edge before adding this one.
	 */
	public boolean constructFromEdge(StsSectionEdge edge)
	{
		try
		{
			if (!addSectionEdge(edge))
			{
				return false;
			}
			int nEdgePoints = edge.getEdgePointsList().getSize();
			boolean isCurved = nEdgePoints > 2;
            return completeSection();
        /*
            geometry = new StsSectionGeometry(this, isCurved, currentModel);
			geometry.computeDipAngle(this);
//			constructReferenceEdge();
			constructSection();
			return true;
	    */
//			constructIntersectedEdges();
		}
		catch (Exception e)
		{
			StsException.outputException("Failed to build section from edge.", e, StsException.WARNING);
			return false;
		}
	}

	public boolean addSectionEdge(StsSectionEdge edge)
	{
		try
		{
			if (edge == null)
			{
				return false;
			}
			byte edgeMatch = checkEdgeMatch(edge);
			if (edgeMatch == MATCH_NOT)
			{
				return false;
			}
			if (edgeMatch == MATCH_REVERSED)
			{
				edge.reverseEdge();

			}
			if (sectionEdges == null)
			{
				sectionEdges = StsObjectRefList.constructor(1, 1, "sectionEdges", this);
			}
			else
			{
				StsXYSurfaceGridable edgeSurface = edge.getSurface();
				int nSectionEdges = sectionEdges.getSize();
				for (int n = nSectionEdges - 1; n >= 0; n--)
				{
					StsSectionEdge existingEdge = (StsSectionEdge) sectionEdges.getElement(n);
					if (existingEdge.getSurface() == edgeSurface)
					{
						existingEdge.delete();
						sectionEdges.delete(n);
					}
				}
			}
			sectionEdges.add(edge);
			edge.setSection(this);
			return true;
		}
		catch (Exception e)
		{
			StsException.outputException("StsSection.addSectionEdge() failed.",
										 e, StsException.WARNING);
			return false;
		}
	}

	private byte checkEdgeMatch(StsSectionEdge edge)
	{
		if (edge == null)
		{
			return MATCH_UNKNOWN;
		}

		StsLine[] lines = edge.getLines();
		if (lines[0] == null || lines[1] == null)
		{
			return MATCH_UNKNOWN;
		}
		StsLine[] sectionLines = getLineList();
		if (sectionLines[0] == null || sectionLines[1] == null)
		{
			return MATCH_UNKNOWN;
		}

		if (lines[0] == sectionLines[0] && lines[1] == sectionLines[1])
		{
			return MATCH;
		}
		if (lines[0] == sectionLines[1] && lines[1] == sectionLines[0])
		{
			return MATCH_REVERSED;
		}
		return MATCH_NOT;
	}

	private StsSurfaceVertex[] getFaultVertices(boolean top)
	{
		if (lines == null || lines.getSize() < 2)
		{
			return null;
		}

		StsSurfaceVertex[] faultVertices = new StsSurfaceVertex[2];
		try
		{
			StsObjectRefList vertices = getFirstLine().getLineVertices();
			faultVertices[0] = top ? (StsSurfaceVertex) vertices.getFirst()
				: (StsSurfaceVertex) vertices.getLast();
			vertices = getLastLine().getLineVertices();
			faultVertices[1] = top ? (StsSurfaceVertex) vertices.getFirst()
				: (StsSurfaceVertex) vertices.getLast();
		}
		catch (NullPointerException e)
		{
			return null;
		}
		return faultVertices;
	}
    /** Each sectionEdge on the section may have a different number of picked points.  Construct a column line
     *  on the section grid for each point on each sectionEdge.  Since all edges have a common first and last column,
     *  the number of columns is 2 + Sum(edges[n].nEdgePoints-2).
     * @return
     */
    private boolean constructCurvedSectionEdges()
	{
        try
        {
            int nEdges = sectionEdges.getSize();
            StsPoint[][] edgePoints = new StsPoint[nEdges][];
            nCols = 2;
            for (int n = 0; n < nEdges; n++)
            {
                StsSectionEdge edge = (StsSectionEdge) sectionEdges.getElement(n);
                edgePoints[n] = edge.getPointsFromEdgePoints();
                int nPoints = edgePoints[n].length;
                nCols += nPoints - 2;
            }

            // edgeGridPoints are edgePoints at each section grid column
            return computeEdgeGridPoints(edgePoints);
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "constructCurvedSectionEdges", e);
            return false;
        }
    }

    public StsPoint[] getSectionIntersectionAtZ(float zCoor)
    {
    	if((zCoor < sectionZMin) || (zCoor > sectionZMax)) return null;
        StsPoint[][] sectionPoints = getPatch().getPoints();
        int row = Math.round((zCoor-sectionZMin)/sectionZInc);
        return sectionPoints[row];
    }

    private boolean constructCurvedSectionGrid()
	{
        try
        {
            StsLine sectionLines[] = getLineList();
            initializePatch(sectionLines); // construct points at first and last columns of grid
            StsPoint[][] sectionPoints = getPatch().getPoints();
            int nEdges = sectionEdges.getSize();
            // edgeGridPoints are edgePoints at each section grid column
            StsPoint[][] edgeGridPoints = new StsPoint[nEdges][];
            for (int e = 0; e < nEdges; e++)
            {
                StsSectionEdge sectionEdge = (StsSectionEdge) sectionEdges.getElement(e);
                edgeGridPoints[e] = sectionEdge.getPoints();
                for (int col = 0; col < nCols; col++)
                {
                    StsPoint edgeGridPoint = edgeGridPoints[e][col];
                    edgeGridPoint.setM(edgeGridPoint.getZorT(isDepth));
                }
            }
            StsLineCoordinates[][] edgePointCoordinates = getLineCoordinates(edgeGridPoints);
            StsLineCoordinates[] rowCoordinates = new StsLineCoordinates[nRows];
            for (int row = 0; row < nRows; row++)
            {
                rowCoordinates[row] = new StsLineCoordinates(sectionPoints[row][0], sectionPoints[row][nCols - 1]);
            }

            for (int col = 1; col < nCols - 1; col++)
            {
                StsPoint[] lineCoordinatePoints = new StsPoint[nEdges];
                for (int e = 0; e < nEdges; e++)
                {
                    lineCoordinatePoints[e] = edgePointCoordinates[e][col].f;
                }
                StsPoint[] lineCoordinateSlopes = StsBezier.computeXYZLineSlopes(lineCoordinatePoints);
                StsPoint[] lineCoordinateColPoints = StsBezier.computeXYZPoints(lineCoordinatePoints, lineCoordinateSlopes, this.sectionZMin, this.sectionZInc, nRows);
                for (int row = 0; row < nRows; row++)
                {
                    sectionPoints[row][col] = rowCoordinates[row].computeTranslatedPoint(lineCoordinateColPoints[row]);
                }
            }
 /*
            StsSeismicVelocityModel velocityModel = currentModel.project.velocityModel;
            if(velocityModel != null)
            {
                for (int col = 0; col < nCols; col++)
                    for (int row = 0; row < nRows; row++)
                        velocityModel.adjustTimeOrDepthPoint(sectionPoints[row][col], isOriginalDepth);
            }
 */
            sectionPatch.makeNormals();
    //        constructCurvedSectionEdges(edgeGridPoints);
            geometryChanged = true;
            return true;
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "constructCurvedSectionGrid", e);
            return false;
        }
    }

    private boolean constructRibbedSectionGrid()
	{
		try
		{
            nCols = lines.getSize();
            sectionPatch = new StsPatch(nRows, nCols, stsColor, this);
			sectionPatch.initialize(sectionZMin, sectionZMax, sectionZInc);
			StsPoint[][] sectionPoints = sectionPatch.getPoints();
			StsPoint point;
            float df = 1.0f/(nCols-1);
            float f = 0.0f;
            for(int col = 0; col < nCols; col++, f += df)
            {
                float z = sectionZMin;
                StsLine line = (StsLine)lines.getElement(col);
                for (int row = 0; row < nRows; row++)
                {
                    point = line.getXYZPointAtZorT(z, true, isDepth);
                    if(point == null) continue;
                    sectionPoints[row][col] = point;
                    point.setF(f);
                    z += sectionZInc;
                }
            }
            /*
            for(int col = 1; col < nCols-1; col++)
            {
                StsLine line = (StsLine)lines.getElement(col);
                if(line.getOnSection() == null)
                    line.setOnSection(this);
            }
            */
            sectionPatch.makeNormals();
		    geometryChanged = true;
            return true;
		}
		catch (Exception e)
		{
			StsException.outputException("StsSection.initializePatch() failed.",
										 e, StsException.WARNING);
			return false;
		}
	}

	public void constructRibbedSectionEdges()
	{
		try
        {
            StsModelSurface[] surfaces = (StsModelSurface[]) currentModel.getCastObjectList(StsModelSurface.class);
            int nSurfaces = surfaces.length;
            if (nSurfaces == 0)
            {
                return;
            }

            if (sectionEdges == null)
            {
                sectionEdges = StsObjectRefList.constructor(nSurfaces, 1, "sectionEdges", this);
            }
            for (int n = 0; n < nSurfaces; n++)
            {
                if (hasSurface(surfaces[n])) continue;
                StsSectionEdge edge = new StsSectionEdge(type, this, surfaces[n], 1);
                sectionEdges.add(edge);
            }
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "constructRibbonSectionEdges", e);
        }
    }
/*
	private void constructCurvedSectionEdges(StsPoint[][] edgeGridPoints)
	{
		if (sectionEdges == null)
		{
			return;
		}
		int nEdges = sectionEdges.getSize();
		for (int n = 0; n < nEdges; n++)
		{
			StsSectionEdge edge = (StsSectionEdge) sectionEdges.getElement(n);
			setEdgeEndCols(edge);
			edge.replaceIntermediateEdgePoints(edgeGridPoints[n]);
			edge.addGridCrossingEdgePoints();
		}
	}
*/
	private void addGridCrossingEdgePoints()
	{
		if (sectionEdges == null) return;
        boolean adjustToGrid = isRibbon() && isVertical();
		int nSectionEdges = sectionEdges.getSize();
        for (int n = 0; n < nSectionEdges; n++)
		{
			StsSectionEdge edge = (StsSectionEdge) sectionEdges.getElement(n);
			edge.addGridCrossingEdgePoints(adjustToGrid);
		}
    }
/*
    private void initializeSectionEdgesTorZ()
    {
		if (sectionEdges == null) return;
		int nSectionEdges = sectionEdges.getSize();
        for (int n = 0; n < nSectionEdges; n++)
		{
			StsSectionEdge edge = (StsSectionEdge) sectionEdges.getElement(n);
            edge.initializeTorZ();
		}
    }
*/
    private boolean initializeSectionGrid()
	{
		if (isRibbon())
		{
			return constructRibbonSectionGrid();
		}
		else if(isCurved())
		{
			return constructCurvedSectionGrid();
		}
        else if(isRibbed())
        {
            return constructRibbedSectionGrid();
        }
        else
            return false;
    }

	private boolean constructRibbonSectionGrid()
	{
        nCols = 2;
        if (!initializePatch(getLineList()))
		{
			return false; // construct points at first and last columns of grid
		}
		sectionPatch.makeNormals();
		geometryChanged = true;
		return true;
	}
/*
	private boolean initializeCurvedSectionGrid()
	{
		try
		{
			if (sectionEdges == null)
			{
				return false;
			}

			StsLine sectionLines[] = getLineList();
			initializePatch(sectionLines); // construct points at first and last columns of grid
			StsPoint[][] sectionPoints = getPatch().getPoints();
			int nEdges = sectionEdges.getSize();
			StsPoint[][] edgeGridPoints = new StsPoint[nEdges][nCols];
			for (int n = 0; n < nEdges; n++)
			{
				StsSectionEdge edge = (StsSectionEdge) sectionEdges.getElement(n);
				StsObjectList edgePoints = edge.getPersistentEdgePoints();
				int nEdgePoints = StsObjectList.getSize(edgePoints);
				if (nEdgePoints != nCols)
				{
					StsException.systemError("StsSection.initializeCurvedSectionGrid() failed." +
											 " nEdgePoints: " + nEdgePoints + " != nSectionCols: " + nCols +
											 " for section: " + getLabel());
					return false;
				}
				for (int col = 0; col < nCols; col++)
				{
					StsGridSectionPoint edgePoint = (StsGridSectionPoint) edgePoints.getElement(col);
					edgeGridPoints[n][col] = edgePoint.getPointXYZorT();
				}
			}

			StsLineCoordinates[][] edgePointCoordinates = getLineCoordinates(edgeGridPoints);
			StsLineCoordinates[] rowCoordinates = new StsLineCoordinates[nRows];
			for (int row = 0; row < nRows; row++)
			{
				rowCoordinates[row] = new StsLineCoordinates(sectionPoints[row][0], sectionPoints[row][nCols - 1]);
			}

			for (int col = 1; col < nCols - 1; col++)
			{
				StsPoint[] lineCoordinatePoints = new StsPoint[nEdges];
				for (int e = 0; e < nEdges; e++)
				{
					lineCoordinatePoints[e] = edgePointCoordinates[e][col].f;
				}
				StsPoint[] lineCoordinateSlopes = StsBezier.computeXYZLineSlopes(lineCoordinatePoints);
				StsPoint[] lineCoordinateColPoints = StsBezier.computeXYZPoints(lineCoordinatePoints, lineCoordinateSlopes, this.sectionZMin, this.dZRow, nRows);
				for (int row = 0; row < nRows; row++)
				{
					sectionPoints[row][col] = rowCoordinates[row].computeTranslatedPoint(lineCoordinateColPoints[row]);
				}
			}

			sectionPatch.makeNormals();
			geometryChanged = true;
			return true;
		}
		catch (Exception e)
		{
			StsException.outputException("StsSection.initializeCurvedSectionGrid() failed.",
										 e, StsException.WARNING);
			return false;
		}
	}
*/
    /** create a section grid column for each interior edge point.
     *  Add a point on each edge at this common column (except for edge which already has that column location).
     *  For each sectionEdge, save the resulting points in the edge.points array.
     */
    private boolean computeEdgeGridPoints(StsPoint[][] edgePoints)
	{
		int e = -1, col = -1;
		int[] indices = null;
		try
		{
			// put arc length in "f" value of each edgePoint
			int nEdges = edgePoints.length;
			for (int n = 0; n < nEdges; n++)
			{
				StsPoint.computeNormalizedArcLengths(StsPoint.XY, edgePoints[n]);
			}
			StsPoint[][] edgeGridPoints = new StsPoint[nEdges][nCols];

			indices = new int[nEdges];
			for (e = 0; e < nEdges; e++)
			{
				edgeGridPoints[e][0] = edgePoints[e][0].copy();
				indices[e] = 1;
			}
			for (col = 1; col < nCols - 1; col++)
			{
				int nextEdgeIndex = -1;
				float minLength = largeFloat;
				for (e = 0; e < nEdges; e++)
				{
					float edgeLength = edgePoints[e][indices[e]].getF();
					if (edgeLength < minLength)
					{
						nextEdgeIndex = e;
						minLength = edgeLength;
					}
				}
				for (e = 0; e < nEdges; e++)
				{
					int index = indices[e];
					if (e == nextEdgeIndex)
					{
						edgeGridPoints[e][col] = edgePoints[e][index].copy();
						indices[e]++;
					}
					else
					{
						StsPoint point0 = edgePoints[e][index - 1];
						StsPoint point1 = edgePoints[e][index];
						float length0 = point0.getF();
						float length1 = point1.getF();
						float f = (minLength - length0) / (length1 - length0);
						edgeGridPoints[e][col] = StsPoint.staticInterpolatePoints(point0, point1, f);
					}
				}
			}
			for (e = 0; e < nEdges; e++)
			{
				int nEdgePoints = edgePoints[e].length;
				edgeGridPoints[e][nCols - 1] = edgePoints[e][nEdgePoints - 1].copy();
			}

            for (e = 0; e < nEdges; e++)
            {
                StsSectionEdge edge = (StsSectionEdge) sectionEdges.getElement(e);
                setEdgeEndCols(edge);
                edge.setPoints(edgeGridPoints[e]);
            }
            return true;
		}
		catch (Exception ex)
		{
			StsException.outputException("StsSection.getEdgeGridPoints() failed. " +
										 " e " + e + " col " + col + " indices[e] " + indices[e] + "edge length " + edgePoints[e].length,
										 ex, StsException.WARNING);
			return false;
		}
	}

	private StsPoint[][] removeCloseEdgeGridPoints(StsPoint[][] edgeGridPoints)
	{
		float ds, dsMin = 0.1f;

		try
		{
			int nRows = edgeGridPoints.length;
			int nCols = edgeGridPoints[0].length;
			int nNewCols = nCols;
			boolean[] remove = new boolean[nCols]; // assume default classInitialize to false
			int prevCol = 0;
			for (int col = 1; col < nCols - 1; col++)
			{
				float s = edgeGridPoints[0][col].getF();
				float sPrev = edgeGridPoints[0][prevCol].getF();
				float sNext = edgeGridPoints[0][col + 1].getF();
				if (sNext == sPrev)
				{
					continue;
				}
				ds = Math.min(s - sPrev, sNext - s) / (sNext - sPrev);
				if (ds <= dsMin)
				{
					remove[col] = true;
					nNewCols--;
				}
				else
				{
					prevCol = col;
				}
			}
			if (nNewCols == nCols)
			{
				return edgeGridPoints;
			}

			StsPoint[][] newEdgeGridPoints = new StsPoint[nRows][nNewCols];
			int nNewCol = 0;
			for (int col = 0; col < nCols; col++)
			{
				if (!remove[col])
				{
					for (int row = 0; row < nRows; row++)
					{
						newEdgeGridPoints[row][nNewCol] = edgeGridPoints[row][col];

					}
					nNewCol++;
				}
			}
			return newEdgeGridPoints;
		}
		catch (Exception e)
		{
			StsException.outputException("StsSection.removeCloseEdgeGridPoints() failed.",
										 e, StsException.WARNING);
			return null;
		}
	}

	private boolean initializePatch(StsLine[] sectionLines)
	{
		try
		{
			sectionPatch = new StsPatch(nRows, nCols, stsColor, this);
			sectionPatch.initialize(sectionZMin, sectionZMax, sectionZInc);
			StsPoint[][] sectionPoints = sectionPatch.getPoints();
			float z = sectionZMin;
			StsPoint point;
			for (int row = 0; row < nRows; row++)
			{
				point = sectionLines[0].getXYZPointAtZorT(z, true, isDepth);
				sectionPoints[row][0] = point;
				point.setF(0.0f);
				point = sectionLines[1].getXYZPointAtZorT(z, true, isDepth);
				sectionPoints[row][nCols - 1] = point;
				point.setF(1.0f);
				z += sectionZInc;
			}
//			sectionPatch.makeNormals();
			return true;
		}
		catch (Exception e)
		{
			StsException.outputException("StsSection.initializePatch() failed.",
										 e, StsException.WARNING);
			return false;
		}
	}

	private StsLineCoordinates[][] getLineCoordinates(StsPoint[][] edgeGridPoints, boolean isDepth)
	{
		int nEdges = edgeGridPoints.length;
		int nCols = edgeGridPoints[0].length;
		StsLineCoordinates[][] lineCoors = new StsLineCoordinates[nEdges][nCols];
		StsLine[] sectionLines = getLineList();
		for (int e = 0; e < nEdges; e++)
		{
			for (int col = 0; col < nCols; col++)
			{
				StsPoint gridPoint = edgeGridPoints[e][col];
				float z = gridPoint.getZorT(isDepth);
				StsPoint point0 = sectionLines[0].getXYZPointAtZorT(z, true, isDepth);
				StsPoint point1 = sectionLines[1].getXYZPointAtZorT(z, true, isDepth);
				lineCoors[e][col] = new StsLineCoordinates(point0, point1, gridPoint, z);
			}
		}
		return lineCoors;
	}

    private StsLineCoordinates[][] getLineCoordinates(StsPoint[][] edgeGridPoints)
	{
		int nEdges = edgeGridPoints.length;
		int nCols = edgeGridPoints[0].length;
		StsLineCoordinates[][] lineCoors = new StsLineCoordinates[nEdges][nCols];
		StsLine[] sectionLines = getLineList();
		for (int e = 0; e < nEdges; e++)
		{
			for (int col = 0; col < nCols; col++)
			{
				StsPoint gridPoint = edgeGridPoints[e][col];
				float z = gridPoint.getZ();
				StsPoint point0 = sectionLines[0].getXYZPointAtZorT(z, true);
				StsPoint point1 = sectionLines[1].getXYZPointAtZorT(z, true);
				lineCoors[e][col] = new StsLineCoordinates(point0, point1, gridPoint, z);
			}
		}
		return lineCoors;
	}

    private StsPoint computeInterpolatedGridPoint(float zRow, int col, StsPoint[][] edgeGridPoints, StsLineCoordinates[][] edgePointCoordinates, StsLineCoordinates newCoor)
	{
		try
		{
			int nEdges = edgeGridPoints.length;

			StsPoint[] edgeGridColPoints = new StsPoint[nEdges];
			for (int e = 0; e < nEdges; e++)
			{
				edgeGridColPoints[e] = edgeGridPoints[e][col];

			}
			float[] weights = getEdgeColWeights(zRow, edgeGridColPoints);

			StsPoint[] translatedPoints = new StsPoint[nEdges];
			StsPoint point = new StsPoint(3);
			point.setZorT(zRow);
			for (int e = 0; e < nEdges; e++)
			{
				StsPoint translatedPoint = StsLineCoordinates.computeTranslatedPoint(edgeGridColPoints[e], edgePointCoordinates[e][col], newCoor);
				point.v[0] += weights[e] * translatedPoint.v[0];
				point.v[1] += weights[e] * translatedPoint.v[1];
//                point.v[2] += weights[e] * translatedPoint.v[2];
//                point.v[4] += weights[e] * translatedPoint.v[4];
			}
			return point;
		}
		catch (Exception e)
		{
			StsException.outputException("StsSection.computeInterpolatedGridPoint() failed.",
										 e, StsException.WARNING);
			return null;
		}
	}

	private float[] getEdgeColWeights(float z, StsPoint[] points)
	{
		int nPoints = points.length;
		float[] weights = new float[nPoints];
		if (nPoints == 1)
		{
			weights[0] = 1.0f;
		}
		else if (nPoints == 2)
		{
			float[] zValues = new float[2];
			for (int n = 0; n < 2; n++)
			{
				zValues[n] = points[n].getZorT();

			}
			weights[1] = (z - zValues[0]) / (zValues[1] - zValues[0]);
			weights[0] = 1.0f - weights[1];
		}
		else
		{
			for (int n = 0; n < nPoints; n++)
			{
				weights[n] = 1.0f / nPoints;

			}
		}
		return weights;
	}

	/** Get tangent on section at point nearest this point.
	 *  @param point point near (or on) section
	 *  @return tangent in horizontal +section direction
	 */
	public StsPoint getTangentAtPoint(StsPoint point)
	{
		return getTangentAtPoint(point, RIGHT);
	}

	public StsPoint getTangentAtPoint(StsPoint point, int side)
	{
		StsPoint tangent;

		if (isRibbon())
		{
			StsPoint[] points = getSectionSidePointsAtZ(point.getZorT(), true);
			tangent = new StsPoint(points[1]);
			tangent.subtract(points[0]);
			tangent.normalize();
		}
		else
		{
			tangent = StsSectionPoint.getTangentAtPoint(point, this);

		}
		if (side == LEFT)
		{
			tangent.reverse();

		}
		return tangent;
	}

	/** This is a general routine used to compute lots of useful things about the
	 *  the relationship of a point to a nearby section.
	 *  A sectionPoint is used to hold all this information:
	 *
	 *  Note that the input value is point and the output value is nearestPoint so
	 *  use them accordingly!!
	 */

	public boolean computeNearestPoint(StsSectionPoint sectionPoint)
//     public StsPoint nearestPoint(StsPoint point, int coorType, float startF, StsPoint tangent,
//                                 StsPoint normal, int[] side, float[] sectionF)
	{
//        if (referenceEdge == null) return false;

		// Since we have only one Section edge, we translated all points to this
		// single coordinate system, compute nearest point and translate back.

		// If we want the nearest point, use getLineIndexNearestPoint.
		// If we want it on a grid X or Y line, use getLineIndexIntersectingLine.
		int gridType = sectionPoint.gridType;

		if (gridType == StsPoint.X || gridType == StsPoint.Y)
		{
			return getLineIndexIntersectingLine(sectionPoint);
		}
		else
		{
			return getLineIndexNearestPoint(sectionPoint);
		}
	}

	public boolean computeNearestPoint(StsSectionPoint sectionPoint, int rowOrCol)
	{
		if (rowOrCol == ROW)
		{
			sectionPoint.gridType = StsPoint.X;
		}
		else if (rowOrCol == COL)
		{
			sectionPoint.gridType = StsPoint.Y;

		}
		return computeNearestPoint(sectionPoint);
	}

	/** Search the line for the point nearest the given point.  Distance is
	 *  computed using the coordinate geomType specified:
	 *      xy  distance is computed in x and y
	 *      xyz distance is computed in x, y, and z
	 *      x   interval must include x, but distance is xy
	 *      y   interval must include y, but distance is xy
	 *  if a reasonable startF is given, this is used to compute a starting
	 *  index and the search is down and up from there
	 */
	private boolean getLineIndexNearestPoint(StsSectionPoint sectionPoint)

	{
		StsPoint[] linePoints = null;
		StsPoint linePoint, linePoint1, linePoint2;
		StsPoint point;
		float distSq1, distSq2, distSq12;
		int minIndex = -1;
		float f;
		float minDistSq;

		if (sectionPatch == null)
		{
			return false;
		}

		point = sectionPoint.point;
		linePoints = sectionPatch.getPointsAtZ(point.getZorT());

		int distanceType = sectionPoint.distanceType;
		float startF = sectionPoint.sectionColF;

		int nLinePoints = linePoints.length;

		if (nLinePoints == 2)
		{
			minIndex = 0;
		}
		else if (startF < -0.0f || startF > nLinePoints - 1) /** startF not used: search whole line */
		{
			minDistSq = largeFloat;

			for (int i = 0; i < nLinePoints; i++)
			{
				float distSq = point.distanceSquaredType(distanceType, linePoints[i]);

				if (distSq < minDistSq)
				{
					minDistSq = distSq;
					minIndex = i;
				}
			}
		}
		else
		/* Use startF to compute a starting index */
		{
			int startI = (int) startF;
			startI = StsMath.minMax(startI, 0, nLinePoints - 1);
			minDistSq = largeFloat;

			for (int i = startI; i < nLinePoints; i++)
			{
				float distSq = point.distanceSquaredType(distanceType, linePoints[i]);

				if (distSq < minDistSq)
				{
					minDistSq = distSq;
					minIndex = i;
				}
				else if (distSq > 2.0f * minDistSq)
				{
					break;
				}
			}

			for (int i = startI; i >= 0; i--)
			{
				float distSq = point.distanceSquaredType(distanceType, linePoints[i]);

				if (distSq < minDistSq)
				{
					minDistSq = distSq;
					minIndex = i;
				}
				else if (distSq > 2.0f * minDistSq)
				{
					break;
				}
			}
		}

		if (minIndex == -1)
		{
			StsException.outputException(new StsException(StsException.WARNING,
				"StsSection.getLineIndexNearestPoint(...) failed.",
				"Section: " + getIndex() + sectionPoint.toString()));
			return false;
		}

		if (minIndex == nLinePoints - 1)
		{
			minIndex--;

		}
		f = computeNearestPoint(point, linePoints, minIndex, distanceType);

		if (f < 0.0f)
		{
			if (minIndex == 0)
			{
				f = 0.0f;
			}
			else
			{
				minIndex--;
				f = computeNearestPoint(point, linePoints, minIndex, distanceType);
				if (f < 0.0f || f > 1.0f)
				{
					minIndex++;
					f = 0.0f;
				}
			}
		}
		else if (f > 1.0f)
		{
			if (minIndex == nLinePoints - 2)
			{
				f = 1.0f;
			}
			else
			{
				minIndex++;
				f = computeNearestPoint(point, linePoints, minIndex, distanceType);
				if (f < 0.0f || f > 1.0f)
				{
					minIndex--;
					f = 1.0f;
				}
			}
		}

		sectionPoint.computeSectionValues(this, minIndex, f, linePoints);
		return true;
	}

	private float computeNearestPoint(StsPoint point, StsPoint[] linePoints, int minIndex, int distanceType)
	{
		StsPoint linePoint1 = linePoints[minIndex];
		StsPoint linePoint2 = linePoints[minIndex + 1];
		float distSq1 = point.distanceSquaredType(distanceType, linePoint1);
		float distSq2 = point.distanceSquaredType(distanceType, linePoint2);
		float distSq12 = linePoint1.distanceSquaredType(distanceType, linePoint2);
		return 0.5f + 0.5f * (distSq1 - distSq2) / distSq12;
	}

	private boolean getLineIndexIntersectingLine(StsSectionPoint sectionPoint)
	{
		int coor;
		float distSq;
		int minIndex = 0;
		float minF = 0.0f;
		float f;
		boolean isIntersected = false;
		StsPoint[] linePoints = null;
		StsPoint linePoint1, linePoint2;

		StsPoint point = sectionPoint.point;

		/** nearestPoint here is a scratch point whose distance from point we minimize */
		StsPoint nearestPoint = new StsPoint(5);

		int gridType = sectionPoint.gridType; /** Indicates which coordinates are variable: X, Y, or XY */
		int distanceType = sectionPoint.distanceType; /** How to measure distance to section */
		float startF = sectionPoint.sectionColF;

		/** EdgePoints is a private class which holds edge geometry at a specific z
		 *  in array linePoints.  If no starting sectionF is available, geometry of
		 *  complete edge is computed; otherwise, points are computed on demand.
		 */

		linePoints = sectionPatch.getPointsAtZ(point.getZorT());
		int nLinePoints = linePoints.length;

		float minDistSq = largeFloat;

		/** If gridType is X, point.y must be in interval (coor = 1). */
		if (gridType == StsPoint.X)
		{
			coor = 1;

			/** If gridType is Y, point.x must be in interval (coor = 0). */
		}
		else if (gridType == StsPoint.Y)
		{
			coor = 0;

			/** If line is in Z-direction, we need to do something different: not
			 *  handled here yet.
			 */
		}
		else
		{
			return false;
		}

		float value = point.v[coor];

		if (startF <= 0.0f || startF > nLinePoints - 1) /** startF not used: search whole line */
		{
			for (int i = 0; i < nLinePoints - 1; i++)
			{
				f = (value - linePoints[i].v[coor]) / (linePoints[i + 1].v[coor] - linePoints[i].v[coor]);

				if (f >= 0.0f && f <= 1.0f)
				{
					nearestPoint.interpolatePoints(linePoints[i], linePoints[i + 1], f);
					distSq = point.distanceSquaredType(distanceType, nearestPoint);

					if (distSq < minDistSq)
					{
						minIndex = i;
						minDistSq = distSq;
						minF = f;
						isIntersected = true;
					}
				}
			}
		}
		else
		/* Use startF to compute a starting index */
		{
			int startI = (int) startF;
			startI = StsMath.minMax(startI, 0, nLinePoints - 2);

			int iUp = startI;
			int iDn = startI - 1;
			boolean iUpOK = iUp >= 0 && iUp < nLinePoints - 1;
			boolean iDnOK = iDn >= 0 && iDn < nLinePoints - 1;

			while (iUpOK || iDnOK)
			{
				if (iUpOK)
				{
					linePoint1 = linePoints[iUp];
					linePoint2 = linePoints[iUp + 1];
					f = (value - linePoint1.v[coor]) / (linePoint2.v[coor] - linePoint1.v[coor]);

					if (f >= 0.0f && f <= 1.0f)
					{
						nearestPoint.interpolatePoints(linePoint1, linePoint2, f);
						distSq = point.distanceSquaredType(distanceType, nearestPoint);

						if (distSq < minDistSq)
						{
							minIndex = iUp;
							minDistSq = distSq;
							minF = f;
							isIntersected = true;
							break;
						}
						else if (distSq > 5.0f * minDistSq)
						{
							iUpOK = false;
						}
					}
				}
				if (iDnOK)
				{
					linePoint1 = linePoints[iDn];
					linePoint2 = linePoints[iDn + 1];
					f = (value - linePoint1.v[coor]) / (linePoint2.v[coor] - linePoint1.v[coor]);

					if (f >= 0.0f && f <= 1.0f)
					{
						nearestPoint.interpolatePoints(linePoint1, linePoint2, f);
						distSq = point.distanceSquaredType(distanceType, nearestPoint);

						if (distSq < minDistSq)
						{
							minIndex = iDn;
							minDistSq = distSq;
							minF = f;
							isIntersected = true;
							break;
						}
						else if (distSq > 5.0f * minDistSq)
						{
							iDnOK = false;
						}
					}
				}
				iUp++;
				iDn--;
				if (iUp < 0 || iUp >= nLinePoints - 1)
				{
					iUpOK = false;
				}
				if (iDn < 0 || iDn >= nLinePoints - 1)
				{
					iDnOK = false;
				}
			}
		}

		if (!isIntersected)
		{
			return false;
		}
		else
		{
			sectionPoint.computeSectionValues(this, minIndex, minF, linePoints);
			return true;
		}
	}

	public boolean setCurrentSection()
	{
		if (this == currentSection)
		{
			return false;
		}
		currentSection = this;
		return true;
	}

	public float getColF(StsPoint point)
	{
		StsSectionPoint sectionPoint = new StsSectionPoint(point);

		if (!computeNearestPoint(sectionPoint))
		{
			return StsParameters.nullValue;
		}
		else
		{
			return sectionPoint.getSectionColF();
		}
	}

	public float getColF(StsPoint point, float guessColF)
	{
		StsSectionPoint sectionPoint = new StsSectionPoint(point);
		sectionPoint.sectionColF = guessColF;

		if (!computeNearestPoint(sectionPoint))
		{
			return StsParameters.nullValue;
		}
		else
		{
			return sectionPoint.getSectionColF();
		}
	}

	/** Project this vertex.point to the nearest point on section; use sectionF as a
	 *  a guess of the the F coordinate.
	 */
	public boolean projectVertexToSection(StsSurfaceVertex vertex)
	{
		StsPoint point = vertex.getPoint();
		StsSectionPoint sectionPoint = new StsSectionPoint(point);
		if (!computeNearestPoint(sectionPoint))
		{
			return false;
		}
		StsGridSectionPoint surfacePoint = vertex.getSurfacePoint();
        surfacePoint.setPoint(sectionPoint.nearestPoint);
        surfacePoint.setRowOrColIndexF(this, COL, sectionPoint.sectionColF);
        vertex.dbFieldChanged("surfacePoint", surfacePoint);
        // surfacePoint.adjustGridRowCol();
        return true;
	}

	public StsPoint getPointOnSection(StsGLPanel3d glPanel3d, StsMousePoint mousePoint)
	{
		if (isRibbon())
		{
			return getPointOnRibbonSection(glPanel3d, mousePoint);
		}
		else if (isCurved())
		{
			return getPointOnCurvedSection(glPanel3d, mousePoint);
		}
		else
		{
			return null;
		}
	}

	public StsPoint getPointOnRibbonSection(StsGLPanel3d glPanel3d, StsMousePoint mousePoint)
	{
		double[] xyz = sectionPatch.getPointOnPatch(glPanel3d, mousePoint);
		StsPoint point = new StsPoint(5);
		point.setX( (float) xyz[0]);
		point.setY( (float) xyz[1]);
		point.setZorT( (float) xyz[2]);

		if (point == null)
		{
			StsMessageFiles.logMessage("Can't find point on section. Try again.");

		}
		return point;
	}

	public StsPoint getPointOnCurvedSection(StsGLPanel3d glPanel3d, StsMousePoint mousePoint)
	{
		StsPoint point = sectionPatch.getPointOnPatchSearch(glPanel3d, mousePoint);

		if (point == null)
		{
			StsMessageFiles.logMessage("Can't find point on section. Try again.");

		}
		return point;
	}

	public StsPoint getEdgeEndPoint(StsSectionEdge sectionEdge, int end)
	{
		// Improve this routine:
		// Instead of getting sectionEdge endPoint which doesn't any section info
		// in it, get an StsSurfaceVertex on the line at that end which has section and f
		if (end == MINUS)
		{
			return sectionEdge.getFirstEdgePoint().getPoint();
		}
		else if (end == PLUS)
		{
			return sectionEdge.getLastEdgePoint().getPoint();
		}
		else
		{
			StsException.outputException(new StsException(StsException.WARNING,
				"StsSection.getEdgeEndPoint(end) failed.", " called with: " + end +
				" value must be MINUS(-1) or PLUS(+1)."));
			return null;
		}
	}

	/** If this line defines the section, return the end it is on or NONE */
	public int getSidePosition(StsLine line)
	{
		if (line == null)
		{
			return NONE;
		}
		if (line == lines.getFirst())
		{
			return MINUS;
		}
		else if (line == lines.getLast())
		{
			return PLUS;
		}
		else
		{
			return NONE;
		}
	}

	public float getColF(StsPoint point, StsLine line)
	{
		if (line != null)
		{
			int linePosition = getSidePosition(line);
			if (linePosition == MINUS)
			{
				return 0.0f;
			}
			else if (linePosition == PLUS)
			{
				return (float) getColMax();
			}
			else
			{
				return getColF(point);
			}
		}
		else
		{
			return getColF(point);
		}
	}

	/** If this line defines the section, return the end it is on as a sectionF
	 *  value: 0.0 for MINUS end and nCols-1 for PLUS end, nullValue otherwise. */
	public float getLineIndexF(StsLine line)
	{
		if (line == lines.getFirst())
		{
			return 0.0f;
		}
		else if (line == lines.getLast())
		{
			return (float) (nCols - 1);
		}
		else
		{
			return StsParameters.nullValue;
		}
	}

	public int getLineIndex(StsLine line)
	{
		if (line == lines.getFirst())
		{
			return 0;
		}
		else if (line == lines.getLast())
		{
			return nCols - 1;
		}
		else
		{
			StsException.systemError("StsSection.getLineIndex() failed." +
									 line.getLabel() + " is not at either end of section: " + getLabel());
			return -1;
		}
	}

	/** Delete this section and the sectionEdge */
	public boolean delete()
	{
		try
		{
			if (sectionEdges != null)
				sectionEdges.forEach("delete");
			if (lines != null)
				lines.forEach("deleteSection", this);
			if (rightLines != null)
				rightLines.forEach("deleteSection", this);
			if (leftLines != null)
				leftLines.forEach("deleteSection", this);

            checkSetModelZDomain();
		}
		catch (Exception e)
		{
			StsException.outputWarningException(this, "delete", e);
			return false;
		}
		try
		{
			super.delete();
			return true;
		}
		catch (Exception e)
		{
			StsException.systemError("StsSection.delete() failed.");
			return false;
		}
	}

    public void deleteSectionEdge(StsModelSurface surface)
    {
        StsObject[] edgeObjects = sectionEdges.getElements();
        for(StsObject edgeObject : edgeObjects)
        {
            StsSectionEdge edge = (StsSectionEdge)edgeObject;
            if(edge.getSurface() == surface)
                sectionEdges.delete(edge);            
        }
    }

    private void checkSetSectionClass()
    {
        if(sectionClass == null) sectionClass = (StsSectionClass)currentModel.getStsClass(StsSection.class);
    }

    private void checkSetModelZDomain()
    {
        checkSetSectionClass();
        sectionClass.checkSetModelZDomain();
    }

    public void deleteFromLines()
	{
		lines.forEach("deleteConnectedSection", this);
	}

	public void deleteLine(StsLine line)
	{
		if (line.getOnSection() != this)
		{
			return;
		}
		int sectionSide = line.getSectionSide();
		if (sectionSide == LEFT)
		{
			leftLines.delete(line);
		}
		else if (sectionSide == RIGHT)
		{
			rightLines.delete(line);
		}
	}

	public StsObjectList addGridCrossingEdgePoints(StsEdge edge, StsObjectList edgePoints, StsXYSurfaceGridable grid)
	{
		StsGridSectionPoint edgePoint0, edgePoint1;
		StsObjectList newEdgePoints;
		StsPoint point;
		StsPoint sectionEdgePoint;
		int extrapolate = NONE;

		try
		{
			edgePoint0 = (StsGridSectionPoint) edgePoints.getFirst();
            edgePoint0 = checkConvertEdgePointToXYZorT(edgePoint0);
            edgePoint1 = (StsGridSectionPoint) edgePoints.getLast();
            edgePoint1 = checkConvertEdgePointToXYZorT(edgePoint1);
            int nEdgePoints = edgePoints.getSize();
			int nEstPoints = StsGridSectionPoint.estimateNCrossingPoints(edgePoint0, edgePoint1, 1.5f, this);
			newEdgePoints = new StsObjectList(nEstPoints, 10);

			for (int n = 1; n < nEdgePoints; n++)
			{
				edgePoint1 = (StsGridSectionPoint) edgePoints.getElement(n);
                edgePoint1 = checkConvertEdgePointToXYZorT(edgePoint1);
                addSectionPointToList(edgePoint0, newEdgePoints);
				StsList gridEdgePoints = StsGridSectionPoint.getGridCrossings(grid, this,
					edgePoint0, edgePoint1, edge, extrapolate, false, MINUS);
				addSectionPointsToList(gridEdgePoints, newEdgePoints);
				edgePoint0 = edgePoint1;
			}
			addSectionPointToList(edgePoint1, newEdgePoints);
			newEdgePoints.trimToSize();
			return newEdgePoints;
		}
		catch (Exception e)
		{
			StsException.outputException("StsSection.constructSurfaceEdgePoints(...sectionEdge) failed.",
										 e, StsException.WARNING);
			return null;
		}
	}

	public StsObjectList constructSurfaceEdgePoints(StsEdge edge, StsXYSurfaceGridable grid, boolean persistent)
	{
		if (isRowColAligned())
		{
			return constructRowColSurfaceEdgePoints(edge, grid, persistent);
		}
		else
		{
			if (isVertical())
			{
				return constructVerticalNonRowColSurfaceEdgePoints(edge, grid, persistent);
			}
			else
			{
				return constructNonRowColSurfaceEdgePoints(edge, grid, persistent);
			}
		}
	}

    private StsGridSectionPoint checkConvertEdgePointToXYZorT(StsGridSectionPoint edgePoint)
    {
        StsPoint point = edgePoint.getPoint();
        if(point.v.length <= 3) return edgePoint;
        StsGridSectionPoint newEdgePoint = edgePoint.copy();
        newEdgePoint.getPoint().convertToXYZorT();
        return newEdgePoint;
    }

    private StsObjectList constructVerticalNonRowColSurfaceEdgePoints(StsEdge edge, StsXYSurfaceGridable grid,
		boolean persistent)
	{
		StsSurfaceVertex vertex0, vertex1;
		StsGridSectionPoint vertexPoint0, vertexPoint1;
		float colf0, colf1, df;
		int nFirst, nLast, inc;
		int nColPoints;
		StsGridSectionPoint colSectionPoint, lastColSectionPoint;
		StsGridPoint gridPoint;
		int extrapolate = NONE;
		int i, n;

		try
		{
			StsXYSurfaceGridable surface = edge.getSurface();
			float minGridSize = getMinGridSize(surface);
			float intervalDistance = 0.5f * minGridSize;

			StsXYSurfaceGridable edgeGrid = (StsXYSurfaceGridable) surface;

			vertex0 = edge.getPrevVertex();
			vertexPoint0 = vertex0.getSurfacePoint();

			vertex1 = edge.getNextVertex();
			vertexPoint1 = vertex1.getSurfacePoint();

			colf0 = vertexPoint0.getColF(this);
			colf1 = vertexPoint1.getColF(this);

			if (colf1 > colf0)
			{
				nFirst = StsMath.above(colf0);
				nLast = StsMath.below(colf1);
//                nFirst = StsMath.ceiling(colf0);
//                nLast = StsMath.floor(colf1);
//                if(nFirst == 0) nFirst++;
//                if(nLast == nCols-1) nLast--;
				inc = 1;
				nColPoints = nLast - nFirst + 1;
			}
			else
			{
				nFirst = StsMath.below(colf0);
				nLast = StsMath.above(colf1);
//                nFirst = StsMath.floor(colf0);
//                nLast = StsMath.ceiling(colf1);
//                if(nLast == 0) nLast++;
//                if(nFirst == nCols-1) nFirst--;
				inc = -1;
				nColPoints = nFirst - nLast + 1;
			}

			int nEstPoints = StsGridSectionPoint.estimateNCrossingPoints(vertexPoint0, vertexPoint1, 1.5f, this);
			nEstPoints += nCols;
			StsObjectList sectionPoints = new StsObjectList(nEstPoints, 10);

			sectionPoints.add(vertexPoint0);

			colSectionPoint = vertexPoint0;

			float startZ = sectionPatch.getZTop();
			for (i = 0, n = nFirst; i < nColPoints; i++, n = n + inc)
			{
				gridPoint = computeSectionColSurfaceIntersect( (float) n, startZ, grid);
				if (gridPoint == null)
				{
					continue;
				}
				lastColSectionPoint = colSectionPoint;
				colSectionPoint = new StsGridSectionPoint(gridPoint, persistent);
				float rowF = getRowF(gridPoint.getZorT());
				float colF = (float) n;
				colSectionPoint = new StsGridSectionPoint(gridPoint, rowF, colF, edge, this, persistent);

// rowColFix		colSectionPoint.setRowOrColIndexF(this, COL, colF);
				getGridCrossings(grid, lastColSectionPoint, colSectionPoint, edge, extrapolate, persistent, MINUS,  sectionPoints);
				sectionPoints.add(colSectionPoint);
//                checkAddSectionPoint(sectionPoints, colSectionPoint);
			}

			lastColSectionPoint = colSectionPoint;
			colSectionPoint = vertexPoint1;
			getGridCrossings(grid, lastColSectionPoint, colSectionPoint, edge, extrapolate, persistent, MINUS, sectionPoints);
			sectionPoints.add(vertexPoint1);

//            checkRemoveClosePoints(sectionPoints);

			sectionPoints.trimToSize();
			return sectionPoints;
		}
		catch (Exception e)
		{
			StsException.outputException("StsSection.edgePointsOnGrid() failed.",
										 e, StsException.WARNING);
			return null;
		}
	}

	// if two points are too close, delete the first and reset the second to the
	// combined row and col position.  If the first is on a row and the second is on a col,
	// for example, set the second to that row.
	private void checkRemoveClosePoints(StsObjectList sectionPoints)
	{
		int nPoints;

		nPoints = sectionPoints.getSize();
		if (nPoints < 3)
		{
			return;
		}

		checkRemoveClosePoints(sectionPoints, nPoints - 1, nPoints - 2);

		nPoints = sectionPoints.getSize();
		if (nPoints < 3)
		{
			return;
		}

		checkRemoveClosePoints(sectionPoints, 0, 1);

		nPoints = sectionPoints.getSize();
		if (nPoints < 3)
		{
			return;
		}

		for (int n = nPoints - 2; n > 0; n--)
		{
			checkRemoveClosePoints(sectionPoints, n - 1, n);
		}
	}

	// check if these two points are too close; if so, adjust the first and delete the second.

	private void checkRemoveClosePoints(StsObjectList sectionPoints, int n0, int n1)
	{
		StsGridSectionPoint point0, point1;
		float r0, c0, r1, c1;
		int rc0, rc1;

		point0 = (StsGridSectionPoint) sectionPoints.getElement(n0);
		r0 = point0.getRowF(null);
		c0 = point0.getColF(null);
		rc0 = point0.getRowOrCol(null);

		point1 = (StsGridSectionPoint) sectionPoints.getElement(n1);
		r1 = point1.getRowF(null);
		c1 = point1.getColF(null);
		rc1 = point1.getRowOrCol(null);

		float dr = r1 - r0;
		float dc = c1 - c0;

		float distSq = dr * dr + dc * dc;
		if (distSq > 2.0e-4)
		{
			return;
		}

		if (rc0 != rc1)
		{
			if ( (rc1 == ROW || rc1 == ROWCOL) && rc0 == COL)
			{
				point0.setRowOrColIndex(null, ROW, (int) r1);
				rc0 = ROWCOL;
				sectionPoints.delete(n1);
				pointDeleted(r0, c0, rc0, r1, c1, rc1, point0);
			}
			else if ( (rc1 == COL || rc1 == ROWCOL) && rc0 == ROW)
			{
				point0.setRowOrColIndex(null, COL, (int) c1);
				rc0 = ROWCOL;
				sectionPoints.delete(n1);
				pointDeleted(r0, c0, rc0, r1, c1, rc1, point0);
			}
		}
		else // same rowOrCol: delete the second
		{
			sectionPoints.delete(n1);
			pointDeleted(r0, c0, rc0, r1, c1, rc1, point0);
		}
	}

	private void pointDeleted(float r0, float c0, int rc0, float r1, float c1, int rc1, StsGridSectionPoint point)
	{
		float newRowF = point.getRowF(null);
		float newColF = point.getColF(null);
		int newRowOrCol = point.getRowOrCol(null);

		StsException.systemDebug("edgePoint deleted: row " + r1 + " col " + c1 + StsParameters.rowCol(rc1) + "\n" +
								 "edgePoint adjusted:  row " + r0 + " col " + c0 + StsParameters.rowCol(rc0) + "\n" +
								 "                to:  row " + newRowF + " col " + newColF +
								 StsParameters.rowCol(newRowOrCol));
	}

	/*
	 private StsObjectList NEWconstructNonRowColSurfaceEdgePoints(StsSectionEdge edge, StsXYSurfaceGridable grid, boolean persistent)
	 {
	  StsSurfaceVertex vertex0, vertex1;
	  StsGridSectionPoint vertexPoint0, vertexPoint1;
	  float sectionColF0, sectionColF1;
	  int nFirst, nLast, inc;
	  float fFirst, fLast;
	  int nColPoints;
	  StsGridSectionPoint colSectionPoint, lastColSectionPoint;
	  StsPoint point;
	  StsGridSectionPoint sectionPoint, lastSectionPoint;
	  StsPoint sectionEdgePoint = new StsPoint();
	  int extrapolate = NONE;
	  int i, n;

	  try
	  {
	   StsGrid surfaceGrid = (StsGrid)(edge.getSurfaceGrid());

	   vertex0 = edge.getPrevVertex();
	   vertexPoint0 = vertex0.getSurfacePoint().getClone();

	   vertex1 = edge.getNextVertex();
	   vertexPoint1 = vertex1.getSurfacePoint().getClone();

	   sectionColF0 = vertexPoint0.getColF(this);
	   sectionColF1 = vertexPoint1.getColF(this);

	   if(sectionColF1 > sectionColF0)
	   {
		nFirst = StsMath.above(sectionColF0);
		nLast = StsMath.below(sectionColF1);
		inc = 1;
		nColPoints = nLast - nFirst + 1;
	   }
	   else
	   {
		nFirst = StsMath.below(sectionColF0);
		nLast = StsMath.above(sectionColF1);
		inc = -1;
		nColPoints = nFirst - nLast + 1;
	   }

	   fFirst = sectionColF0;
	   fLast = sectionColF1;

	   int nEstPoints = StsGridSectionPoint.estimateNCrossingPoints(vertexPoint0, vertexPoint1, 1.5f, this);
	   StsObjectList sectionPoints = new StsObjectList(nEstPoints, 10);

	   sectionPoints.add(vertexPoint0);

	   colSectionPoint = vertexPoint0;
	   sectionColF1 = fFirst;
	   float guessZ = vertexPoint0.getPoint().getZ();
	   for(i = 1, n = nFirst; i <= nColPoints; i++, n = n + inc)
	   {
		point = computeSectionColSurfaceIntersect((float)n, guessZ, grid);
		if(point == null) continue;
		guessZ = point.getZ();
		lastColSectionPoint = colSectionPoint;
		sectionColF0 = sectionColF1;
		sectionColF1 = (float)n;
	 float rowF = getRowF(guessZ);
	 colSectionPoint = new StsGridSectionPoint(point, rowF, sectionColF1, (StsEdgeLinkable)edge, this, persistent);
//rowColFix			    colSectionPoint.setRowOrColIndexF(this, COL, sectionColF1);
	 NEWaddIntermediateSectionPoints(lastColSectionPoint, colSectionPoint, sectionColF0, sectionColF1, edge,
		   grid, sectionPoints, persistent);
		sectionPoints.add(colSectionPoint);
	   }

	   lastColSectionPoint = colSectionPoint;

	   colSectionPoint = vertexPoint1;
	   sectionColF0 = sectionColF1;
	   sectionColF1 = fLast;

	   NEWaddIntermediateSectionPoints(lastColSectionPoint, colSectionPoint, sectionColF0, sectionColF1, edge,
		   grid, sectionPoints, persistent);

	   StsGridSectionPoint lastPoint = (StsGridSectionPoint)sectionPoints.getLast();
	   if(vertexPoint1 != lastPoint) sectionPoints.add(vertexPoint1);

	   return sectionPoints;
	  }
	  catch(Exception e)
	  {
	   StsException.outputException("StsSection.edgePointsOnGrid() failed.",
		e, StsException.WARNING);
	   return null;
	  }
	 }
	 */
	/*
	 private StsPoint XcomputeSectionColSurfaceIntersect(float sectionColF, float startZ, StsXYSurfaceGridable grid, boolean extrapolate)
	 {
	  float sectionRowF;
	  float[] xyz = null;
	  int iter = 0;
	  float surfZ, colZ;
	  float errorZ = largeFloat;

	  try
	  {
	   surfZ = startZ;

	   while(iter++ < maxIter)
	   {
		colZ = surfZ;
		sectionRowF = getRowF(colZ);
		xyz = sectionPatch.getPoint(sectionRowF, sectionColF);
		if(xyz == null) return null;
		surfZ = grid.interpolateBilinearZ(xyz, true, true);
		if(surfZ == nullValue) return null;

		errorZ = Math.abs(colZ - surfZ);
		if(errorZ < maxError) return new StsPoint(xyz);
	   }

	 StsException.systemError("StsSection.computeSectionColSurfaceIntersect() failed to converge for: " + getLabel() +
	 " colF: " + sectionColF + " surface: " + grid.getLabel() + ".\n Iterations: " + iter + " errorZ: " + errorZ);
	   return new StsPoint(xyz);
	  }
	  catch(Exception e)
	  {
	 StsException.outputException("StsSection.computeSectionColSurfaceIntersect() failed for: " + getLabel() +
		" colF: " + sectionColF + " surface: " + grid.getLabel() + ".\n Iterations: " + iter + " errorZ: " + errorZ, e, StsException.WARNING);
	   return null;
	  }
	 }
	 */
	private StsGridPoint computeSectionColSurfaceIntersect(float sectionColF, float startZ, StsXYSurfaceGridable grid)
	{
		int rowAbove, rowBelow;
		boolean aboveOK, belowOK;
		StsGridPoint point;
		float surfZ;

		try
		{
			rowAbove = (int) getRowF(startZ);
			rowAbove = StsMath.minMax(rowAbove, 0, nRows - 2);
			rowBelow = rowAbove + 1;
			aboveOK = true;
			belowOK = rowBelow < nRows - 2;

			while (aboveOK || belowOK)
			{
				if (aboveOK)
				{
					point = computeGridIntervalIntersect(grid, rowAbove, sectionColF);
					if (point != null)
					{
						return point;
					}
					if (rowAbove-- <= 0)
					{
						aboveOK = false;
					}
				}
				if (belowOK)
				{
					point = computeGridIntervalIntersect(grid, rowBelow - 1, sectionColF);
					if (point != null)
					{
						return point;
					}
					if (rowBelow++ == nRows - 2)
					{
						belowOK = false;
					}
				}
			}
			// didn't intersect: use top or bot point, whichever is closer
			StsPoint topPoint = sectionPatch.getStsPoint(0.0f, sectionColF);
            float topZ = topPoint.getZorT();
            StsGridPoint topGridPoint = new StsGridPoint(topPoint, grid);
			surfZ = grid.interpolateBilinearZ(topGridPoint, true, true);
			float topDZ = Math.abs(surfZ - topZ);
			StsPoint botPoint = sectionPatch.getStsPoint( (float) (nRows - 1), sectionColF);
            float botZ = botPoint.getZorT();
            StsGridPoint botGridPoint = new StsGridPoint(botPoint, grid);
			surfZ = grid.interpolateBilinearZ(botGridPoint, true, true);
			float botDZ = Math.abs(surfZ - botZ);
			if (topDZ <= botDZ)
			{
				return topGridPoint;
			}
			else
			{
				return botGridPoint;
			}
		}
		catch (Exception e)
		{
			StsException.outputException("StsSection.computeGridIntersect() failed." +
										 " col: " + sectionColF + getLabel(),
										 e, StsException.WARNING);
			return null;
		}
	}

	private StsGridPoint computeGridIntervalIntersect(StsXYSurfaceGridable grid, int rowAbove, float sectionColF)
	{
		StsPoint point0 = sectionPatch.getStsPoint( (float) rowAbove, sectionColF);
		StsPoint point1 = sectionPatch.getStsPoint( (float) (rowAbove + 1), sectionColF);
		StsGridCrossings gridCrossings = new StsGridCrossings(point0, point1, grid);
		return gridCrossings.getGridIntersection(grid);
	}

	/*
	 private StsPoint computeSectionColSurfaceIntersect(float sectionColF, float startZ, StsXYSurfaceGridable grid, boolean extrapolate)
	 {
	  float sectionRowF;
	  float[] xyz = null;
	  int iter = 0;
	  float surfZ, colZ;
	  float errorZ = largeFloat;

	  try
	  {
	   surfZ = startZ;

	   while(iter++ < maxIter)
	   {
		colZ = surfZ;
		sectionRowF = getRowF(colZ);
		xyz = sectionPatch.getPoint(sectionRowF, sectionColF);
		if(xyz == null) return null;
		surfZ = grid.interpolateBilinearZ(xyz, true, true);
		if(surfZ == nullValue) return null;

		errorZ = Math.abs(colZ - surfZ);
		if(errorZ < maxError) return new StsPoint(xyz);
	   }

	   // failed to find an intersection iteratively: search in opposite direction sequentially
	   float failedDirection = surfZ - startZ;
	   return computeSectionColSurfaceIntersect(sectionColF, startZ, failedDirection, grid, extrapolate);
	  }
	  catch(Exception e)
	  {
	 StsException.outputException("StsSection.computeSectionColSurfaceIntersect() failed for: " + getLabel() +
		" colF: " + sectionColF + " surface: " + grid.getLabel() + ".\n Iterations: " + iter + " errorZ: " + errorZ, e, StsException.WARNING);
	   return null;
	  }
	 }

// search sequentially in opposite direction
	 private StsPoint computeSectionColSurfaceIntersect(float sectionColF, float startZ, float failedDirection, StsXYSurfaceGridable grid, boolean extrapolate)
	 {
	  StsPoint intersection;
	  int searchDirection;

	  if(failedDirection >= 0.0f)
	   searchDirection = -1;
	  else
	   searchDirection = 1;

	  intersection = computeSectionColSurfaceIntersect(sectionColF, startZ, searchDirection, grid, extrapolate);
	  if(intersection != null) return intersection;
	  intersection = computeSectionColSurfaceIntersect(sectionColF, startZ, -searchDirection, grid, extrapolate);
	  if(intersection != null) return intersection;

	  StsException.systemError("StsSection.computeSectionColSurfaceIntersect() failed to converge for: " + getLabel() +
	   " colF: " + sectionColF + " surface: " + grid.getLabel());
	  return null;
	 }

	 private StsPoint computeSectionColSurfaceIntersect(float sectionColF, float startZ, int searchDirection, StsXYSurfaceGridable grid, boolean extrapolate)
	 {
	  float sectionRowF;
	  int startRow, endRow, incRow;
	  int nSearchRows;
	  float[] xyz = null;
	  float surfZ, colZ;
	  float nextDZ, prevDZ;
	  float prevSurfZ;
	  float debugRowF, debugColF;
	  try
	  {
	   sectionRowF = getRowF(startZ);
	   float zInc = currentModel.getProject().getZInc();
	   if(searchDirection >= 0)
	   {
		startRow = StsMath.floor(sectionRowF);
	 endRow = nRows-1;
	 incRow = 1;
	 nSearchRows = endRow - startRow;
	   }
	   else
	   {
	 startRow = StsMath.ceiling(sectionRowF);
	 endRow = 0;
	 incRow = -1;
	 nSearchRows = startRow - endRow;
	 zInc = -zInc;
	   }
	   if(nSearchRows < 1) return null;

	   colZ = getZAtRow(startRow);
	   xyz = sectionPatch.getPoint((float)startRow, sectionColF);
	   if(xyz == null) return null;
//			debugRowF = grid.getRowCoor(xyz[1]);
//			debugColF = grid.getColCoor(xyz[0]);
	   surfZ = grid.interpolateBilinearZ(xyz, true, true);

	   nextDZ = surfZ - colZ;
	   int row = startRow;
	   for(int n = 0; n < nSearchRows; n++)
	   {
	 prevDZ = nextDZ;
	 prevSurfZ = surfZ;
	 row += incRow;
	 colZ += zInc;
		xyz = sectionPatch.getPoint((float)row, sectionColF);
		if(xyz == null) return null;
//			    debugRowF = grid.getRowCoor(xyz[1]);
//			    debugColF = grid.getColCoor(xyz[0]);
		surfZ = grid.interpolateBilinearZ(xyz, true, true);
	 nextDZ = surfZ - colZ;

		if(nextDZ == 0.0f)
		 return new StsPoint(xyz);

	 if(nextDZ > 0.0f && prevDZ <= 0.0f ||
		nextDZ < 0.0f && prevDZ >= 0.0f )
	 {
		 float f = prevDZ/(prevDZ - nextDZ);
	  surfZ = prevSurfZ + f*(surfZ - prevSurfZ);
	  sectionRowF = getRowF(surfZ);
		 xyz = sectionPatch.getPoint(sectionRowF, sectionColF);
//			        debugRowF = grid.getRowCoor(xyz[1]);
//			        debugColF = grid.getColCoor(xyz[0]);
//                    surfZ = grid.interpolateBilinearZ(xyz, true, true);
//				    colZ = xyz[2];
//				    float errorZ = colZ - surfZ;
	  return new StsPoint(xyz);
	 }
	   }
	   return null;
	  }
	  catch(Exception e)
	  {
	 StsException.outputException("StsSection.computeSectionColSurfaceIntersect() failed for: " + getLabel() +
		" colF: " + sectionColF + " surface: " + grid.getLabel(), e, StsException.WARNING);
	   return null;
	  }
	 }
	 */
	private StsObjectList constructNonRowColSurfaceEdgePoints(StsEdge edge, StsXYSurfaceGridable grid, boolean persistent)
	{
		int nFirst, nLast, inc;
		float fFirst, fLast;
		int nColPoints;
		StsGridSectionPoint colSectionPoint, lastColSectionPoint;
		StsGridPoint gridPoint;
		int extrapolate = NONE;
		int i, n;

		try
		{
			StsXYSurfaceGridable surface = edge.getSurface();
			float minGridSize = getMinGridSize(surface);
			float intervalDistance = 0.5f * minGridSize;

			StsSurfaceVertex vertex0 = edge.getPrevVertex();
            StsGridSectionPoint vertexPoint0 = new StsGridSectionPoint(vertex0.getSurfacePoint(), persistent);

			StsSurfaceVertex vertex1 = edge.getNextVertex();
            StsGridSectionPoint vertexPoint1 = new StsGridSectionPoint(vertex1.getSurfacePoint(), persistent);

			float colf0 = vertexPoint0.getColF(this);
			float colf1 = vertexPoint1.getColF(this);

			if (colf1 > colf0)
			{
				nFirst = StsMath.above(colf0);
				nLast = StsMath.below(colf1);
				inc = 1;
				nColPoints = nLast - nFirst + 1;
			}
			else
			{
				nFirst = StsMath.below(colf0);
				nLast = StsMath.above(colf1);
				inc = -1;
				nColPoints = nFirst - nLast + 1;
			}

			fFirst = colf0;
			fLast = colf1;

			int nEstPoints = StsGridSectionPoint.estimateNCrossingPoints(vertexPoint0, vertexPoint1, 1.5f, this);
			StsObjectList sectionPoints = new StsObjectList(nEstPoints, 10);

			sectionPoints.add(vertexPoint0);

			colSectionPoint = vertexPoint0;
			colf1 = fFirst;
			float guessZ = vertexPoint0.getPoint().getZorT();
			for (i = 1, n = nFirst; i <= nColPoints; i++, n = n + inc)
			{
				gridPoint = computeSectionColSurfaceIntersect( (float) n, guessZ, grid);
                if (gridPoint == null)
				{
					continue;
				}
				guessZ = gridPoint.getZorT(isDepth);

				lastColSectionPoint = colSectionPoint;
				colf0 = colf1;
				colf1 = (float) n;
				float rowF = getRowF(guessZ);
				colSectionPoint = new StsGridSectionPoint(gridPoint, rowF, colf1, edge, this, persistent);
// rowColFix        colSectionPoint.setRowOrColIndexF(this, COL, colf1);
				addIntermediateSectionPoints(lastColSectionPoint, colSectionPoint, colf0, colf1, edge,
											 intervalDistance, grid, sectionPoints, extrapolate, persistent);
//                StsList gridEdgePoints = StsGridSectionPoint.getGridCrossings(null, this,
//							    lastColSectionPoint, colSectionPoint, edge, extrapolate, false, MINUS);
//                addSectionPointsToList(gridEdgePoints, sectionPoints);
				sectionPoints.add(colSectionPoint);
			}

			lastColSectionPoint = colSectionPoint;

			colSectionPoint = vertexPoint1;
			colf0 = colf1;
			colf1 = fLast;

			addIntermediateSectionPoints(lastColSectionPoint, colSectionPoint, colf0, colf1, edge,
										 intervalDistance, grid, sectionPoints, extrapolate, persistent);
//            StsList gridEdgePoints = StsGridSectionPoint.getGridCrossings(null, this,
//							lastColSectionPoint, colSectionPoint, edge, extrapolate, false, MINUS);
//            addSectionPointsToList(gridEdgePoints, sectionPoints);

			StsGridSectionPoint lastPoint = (StsGridSectionPoint) sectionPoints.getLast();
			if (vertexPoint1 != lastPoint)
			{
				sectionPoints.add(vertexPoint1);

			}
			sectionPoints = addSectionRowCrossings(sectionPoints, edge);
			sectionPoints.trimToSize();
			return sectionPoints;
		}
		catch (Exception e)
		{
			StsException.outputException("StsSection.edgePointsOnGrid() failed.",
										 e, StsException.WARNING);
			return null;
		}
	}


    public float getMinGridSize(StsXYSurfaceGridable surface)
    {
        return Math.min(surface.getXInc(), surface.getYInc());
    }

    private void addIntermediateSectionPoints(StsGridSectionPoint sectionPoint0,
											  StsGridSectionPoint sectionPoint1, float colf0, float colf1,
											  StsEdgeLinkable edge,
											  float intervalDistance, StsXYSurfaceGridable grid,
											  StsObjectList sectionPoints, int extrapolate, boolean persistent)
	{
		StsGridSectionPoint sectionPoint;

		try
		{
			float distance = sectionPoint0.point.distanceType(StsPoint.XY, sectionPoint1.point);
			int nIntervals = (int) (distance / intervalDistance);
			if (nIntervals > 1)
			{
				float dff = (colf1 - colf0) / nIntervals;
				float colfi = colf0 + dff;
				float guessZ = sectionPoint0.point.getZorT();
				for (int ni = 0; ni < nIntervals - 1; ni++)
				{
					StsGridPoint gridPoint = computeSectionColSurfaceIntersect(colfi, guessZ, grid);
					if (gridPoint != null)
					{
						guessZ = gridPoint.getZorT();
						float rowF = getRowF(guessZ);

						sectionPoint = new StsGridSectionPoint(gridPoint, rowF, colfi, edge, this, false);
						StsList gridEdgePoints = StsGridSectionPoint.getGridCrossings(null, this,
							sectionPoint0, sectionPoint, edge, extrapolate, false, MINUS);
						addSectionPointsToList(gridEdgePoints, sectionPoints);
//                         addRowSectionPoints(sectionPoint0, sectionPoint, edge, grid, sectionPoints, NONE, persistent);
						sectionPoint0 = sectionPoint;
					}
					colfi += dff;
				}
				StsList gridEdgePoints = StsGridSectionPoint.getGridCrossings(null, this,
					sectionPoint0, sectionPoint1, edge, extrapolate, false, MINUS);
				addSectionPointsToList(gridEdgePoints, sectionPoints);
//                addRowSectionPoints(sectionPoint0, sectionPoint1, edge, grid, sectionPoints, extrapolate, persistent);
			}
			else
			{
				StsList gridSectionPoints = StsGridSectionPoint.getGridCrossings(null, this,
					sectionPoint0, sectionPoint1, edge, extrapolate, false, MINUS);
				addSectionPointsToList(gridSectionPoints, sectionPoints);
//               addRowSectionPoints(sectionPoint0, sectionPoint1, edge, grid, sectionPoints, extrapolate, persistent);
			}
		}
		catch (Exception e)
		{
			StsException.outputException("StsSection.addIntermediateColSectionPoints() failed.",
										 e, StsException.WARNING);
		}
	}

	private StsObjectList addSectionRowCrossings(StsObjectList sectionPoints, StsEdgeLinkable edge)
	{
		StsGridSectionPoint sectionPoint0, sectionPoint1;
		int nSectionPoints = sectionPoints.getSize();

		StsObjectList newSectionPoints = new StsObjectList(nSectionPoints, 100);

		sectionPoint1 = (StsGridSectionPoint) sectionPoints.getElement(0);
		for (int n = 1; n < nSectionPoints; n++)
		{
			sectionPoint0 = sectionPoint1;
			sectionPoint1 = (StsGridSectionPoint) sectionPoints.getElement(n);

			addSectionPointToList(sectionPoint0, newSectionPoints);

			StsObjectList gridSectionPoints = getSectionRowCrossings(sectionPoint0, sectionPoint1, edge);
			addSectionPointsToList(gridSectionPoints, newSectionPoints);
		}
		addSectionPointToList(sectionPoint1, newSectionPoints);
		return newSectionPoints;
	}

	private StsObjectList getSectionRowCrossings(StsGridSectionPoint sectionPoint0,
												 StsGridSectionPoint sectionPoint1, StsEdgeLinkable edge)
	{
		float rowF0, rowF1;
		int rowStart, rowEnd, dRow, row;
		int nPoints;
		StsGridSectionPoint sectionPoint;

		try
		{
			rowF0 = sectionPoint0.getRowF(this);
			rowF1 = sectionPoint1.getRowF(this);
         /*
			if (rowF0 < 0 || rowF0 > nRows - 1 || rowF1 < 0 || rowF1 > nRows - 1)
			{
				return null;
			}
         */
			if (rowF1 == rowF0)
			{
				return null;
			}
			else if (rowF1 > rowF0)
			{
				rowStart = StsMath.ceiling(rowF0);
				rowEnd = StsMath.floor(rowF1);
				if (rowStart > rowEnd)
				{
					return null;
				}
				dRow = 1;
				nPoints = rowEnd - rowStart + 1;
			}
			else // rowF1 < rowF0
			{
				rowStart = StsMath.floor(rowF0);
				rowEnd = StsMath.ceiling(rowF1);
				if (rowStart < rowEnd)
				{
					return null;
				}
				dRow = -1;
				nPoints = rowStart - rowEnd + 1;
			}

			StsObjectList newSectionPoints = new StsObjectList(nPoints);

			row = rowStart;
			while (true)
			{
				float ff = (row - rowF0) / (rowF1 - rowF0);
				if (ff > 1.0f)
				{
					break;
				}
				sectionPoint = StsGridSectionPoint.sectionInterpolate(sectionPoint0, sectionPoint1, (double) ff, edge, this, false);
				sectionPoint.setRowOrColIndex(this, ROW, row, false);
				newSectionPoints.add(sectionPoint);
				if (row == rowEnd)
				{
					break;
				}
				row += dRow;
			}
			return newSectionPoints;
		}
		catch (Exception e)
		{
			StsException.outputException("StsSection.getSectionRowCrossings() failed.",
										 e, StsException.WARNING);
			return null;
		}
	}

	/*
	 private void NEWaddIntermediateSectionPoints(StsGridSectionPoint sectionPoint0,
	 StsGridSectionPoint sectionPoint1, float sectionColF0, float sectionColF1, StsEdgeLinkable edge,
		 StsXYSurfaceGridable grid, StsObjectList sectionPoints, boolean persistent)
	 {
	  StsGridSectionPoint sectionPoint;

	  try
	  {
	   sectionPoint = sectionPoint0;
	   while(true)
	   {
	 sectionPoint = getNextGridCrossing(sectionPoint, sectionColF0, sectionColF1, edge, grid, persistent);
		if(sectionPoint == null) return;
		sectionPoints.add(sectionPoint);
	   }
	  }
	  catch(Exception e)
	  {
	   StsException.outputException("StsSection.addIntermediateColSectionPoints() failed.",
		e, StsException.WARNING);
	  }
	 }

	 private StsGridSectionPoint getNextGridCrossing(StsGridSectionPoint sectionPoint,
			 float sectionColF0, float sectionColF1,
			 StsEdgeLinkable edge, StsXYSurfaceGridable grid, boolean persistent)
	 {
	  float lineGridRowF0, lineGridColF0, lineGridRowF1, lineGridColF1;
	  int nextGridRow = -1, nextGridCol = -1;
	  boolean rowIsNext, colIsNext;
	  float gridRowF, gridColF;

	  try
	  {
	   int[] gridRowCol = sectionPoint.getLowerLeftRowCol(null);
	   int gridRow = gridRowCol[0];
	   int gridCol = gridRowCol[1];

	   float z = sectionPoint.getZ();
	   int iter = 0;
	   float errorZ = largeFloat;
	   while(iter++ < maxIter)
	   {
		float sectionRowF = getRowF(z);
		float[] xyz0 = sectionPatch.getPoint(sectionRowF, sectionColF0);
		float[] xyz1 = sectionPatch.getPoint(sectionRowF, sectionColF1);

		lineGridRowF0 = grid.getRowCoor(xyz0[1]);
		lineGridColF0 = grid.getColCoor(xyz0[0]);
		lineGridRowF1 = grid.getRowCoor(xyz1[1]);
		lineGridColF1 = grid.getColCoor(xyz1[0]);

		if(lineGridRowF1 == lineGridRowF0 && lineGridColF1 == lineGridColF0) return null;

		if(lineGridRowF1 > lineGridRowF0)
		{
		 nextGridRow = StsMath.above(gridRow);
		 rowIsNext = nextGridRow < lineGridRowF1;
		}
		else if(lineGridRowF1 < lineGridRowF0)
		{
		 nextGridRow  = StsMath.below(gridRow);
		 rowIsNext = nextGridRow > lineGridRowF1;
		}
		else
		 rowIsNext = false;

		if(lineGridColF1 > lineGridColF0)
		{
		 nextGridCol = StsMath.above(gridCol);
		 colIsNext = nextGridCol < lineGridColF1;
		}
		else if(lineGridColF1 < lineGridColF0)
		{
		 nextGridCol  = StsMath.below(gridCol);
		 colIsNext = nextGridCol > lineGridColF1;
		}
		else
		 colIsNext = false;

		float fRow = largeFloat;
		float fCol = largeFloat;
		float sectionColF;

		if(rowIsNext)
		{
		 rowIsNext = StsMath.between(nextGridRow, lineGridRowF0, lineGridRowF1);
		 if(rowIsNext) fRow = (nextGridRow - lineGridRowF0)/(lineGridRowF1 - lineGridRowF0);
		}
		if(colIsNext)
		{
		 colIsNext = StsMath.between(nextGridCol, lineGridColF0, lineGridColF1);
		 if(colIsNext) fCol = (nextGridCol - lineGridColF0)/(lineGridColF1 - lineGridColF0);
		}

		if(!rowIsNext && !colIsNext)
		 return null;

		float[] xyz;

		if(rowIsNext && colIsNext)
		{
		 if(fRow <= fCol)
		  colIsNext = false;
		 else
		  rowIsNext = false;
		}

		if(rowIsNext)
		 xyz = StsMath.interpolate(xyz0, xyz1, fRow);
		else // colIsNext
		 xyz = StsMath.interpolate(xyz0, xyz1, fCol);

		float gridZ = grid.interpolateBilinearZ(xyz, true, true);
		errorZ = Math.abs(gridZ - z);
		z = gridZ;
		if(errorZ < maxError)
		{
		 if(rowIsNext)
		 {
		  gridColF = lineGridColF0 + fRow*(lineGridColF1 - lineGridColF0);
		  gridRowF = (float)nextGridRow;
		  sectionColF = sectionColF0 + fRow*(sectionColF1 - sectionColF0);
		 }
		 else // colIsNext
		 {
		  gridRowF = lineGridRowF0 + fCol*(lineGridRowF1 - lineGridRowF0);
		  gridColF = (float)nextGridCol;
		  sectionColF = sectionColF0 + fCol*(sectionColF1 - sectionColF0);
		 }

		 StsPoint point = new StsPoint(xyz);
		 sectionPoint = new StsGridSectionPoint(point, gridRowF, gridColF, edge, null, persistent);
		 sectionPoint.addSectionRowCol(sectionRowF, sectionColF, this);
//rowColFix
//					if(rowIsNext)
//						sectionPoint.setRowOrColIndexF(null, ROW, gridRowF);
//					if(colIsNext)
//						sectionPoint.setRowOrColIndexF(null, COL, gridColF);

		 return sectionPoint;
		}
	   }
	   return null;
	  }
	  catch(Exception e)
	  {
	   StsException.outputException("StsSection.getNextGridCrossing() failed.",
		e, StsException.WARNING);
	   return null;
	  }
	 }
	 */
	/*
	 private StsGridSectionPoint XgetNextGridCrossing(StsGridSectionPoint sectionPoint0,StsGridSectionPoint sectionPoint1,
			 float sectionColF0, float sectionColF1,
			 StsEdgeLinkable edge, StsXYSurfaceGridable grid, boolean persistent)
	 {
	  float gridRowF0, gridColF0, gridRowF1, gridColF1;
	  int nextGridRow = -1, nextGridCol = -1;
	  boolean rowIsNext, colIsNext;
	  float gridRowF, gridColF;

	  try
	  {
	   gridRowF0 = sectionPoint0.getRowF(null);
	   gridColF0 = sectionPoint0.getColF(null);
	   gridRowF1 = sectionPoint1.getRowF(null);
	   gridColF1 = sectionPoint1.getColF(null);

	   if(gridRowF1 == gridRowF0 && gridColF1 == gridColF0) return null;

	   if(gridRowF1 > gridRowF0)
	   {
		nextGridRow = StsMath.above(gridRowF0);
		rowIsNext = nextGridRow < gridRowF1;
	   }
	   else if(gridRowF1 < gridRowF0)
	   {
		nextGridRow  = StsMath.below(gridRowF0);
		rowIsNext = nextGridRow > gridRowF1;
	   }
	   else
		rowIsNext = false;

	   if(gridColF1 > gridColF0)
	   {
		nextGridCol = StsMath.above(gridColF0);
		colIsNext = nextGridCol < gridColF1;
	   }
	   else if(gridRowF1 < gridRowF0)
	   {
		nextGridCol  = StsMath.below(gridColF0);
		colIsNext = nextGridCol > gridColF1;
	   }
	   else
		colIsNext = false;

	   float z = sectionPoint0.getZ();
	   int iter = 0;
	   float errorZ = largeFloat;
	   while(iter++ < maxIter)
	   {
		float sectionRowF = getRowF(z);
		float[] xyz0 = sectionPatch.getPoint(sectionRowF, sectionColF0);
		float[] xyz1 = sectionPatch.getPoint(sectionRowF, sectionColF1);

		float lineGridRowF0 = grid.getRowCoor(xyz0[1]);
		float lineGridColF0 = grid.getColCoor(xyz0[0]);
		float lineGridRowF1 = grid.getRowCoor(xyz1[1]);
		float lineGridColF1 = grid.getColCoor(xyz1[0]);

		float fRow = largeFloat;
		float fCol = largeFloat;
		float sectionColF;

		if(rowIsNext)
		{
		 rowIsNext = StsMath.between(nextGridRow, lineGridRowF0, lineGridRowF1);
		 if(rowIsNext) fRow = (nextGridRow - lineGridRowF0)/(lineGridRowF1 - lineGridRowF0);
		}
		if(colIsNext)
		{
		 colIsNext = StsMath.between(nextGridCol, lineGridColF0, lineGridColF1);
		 if(colIsNext) fCol = (nextGridCol - lineGridColF0)/(lineGridColF1 - lineGridColF0);
		}

		if(!rowIsNext && !colIsNext)
		 return null;

		float[] xyz;

		if(rowIsNext && colIsNext)
		{
		 if(fRow <= fCol)
		  colIsNext = false;
		 else
		  rowIsNext = false;
		}

		if(rowIsNext)
		 xyz = StsMath.interpolate(xyz0, xyz1, fRow);
		else // colIsNext
		 xyz = StsMath.interpolate(xyz0, xyz1, fCol);

		float gridZ = grid.interpolateBilinearZ(xyz, true, true);
		errorZ = Math.abs(gridZ - z);
		if(errorZ < maxError)
		{
		 if(rowIsNext)
		 {
		  gridColF = lineGridColF0 + fRow*(lineGridColF1 - lineGridColF0);
		  gridRowF = (float)nextGridRow;
		  sectionColF = sectionColF0 + fRow*(sectionColF1 - sectionColF0);
		 }
		 else // colIsNext
		 {
		  gridRowF = lineGridRowF0 + fCol*(lineGridRowF1 - lineGridRowF0);
		  gridColF = (float)nextGridCol;
		  sectionColF = sectionColF0 + fCol*(sectionColF1 - sectionColF0);
		 }

		 StsPoint point = new StsPoint(xyz);
	 StsGridSectionPoint sectionPoint = new StsGridSectionPoint(point, gridRowF, gridColF, edge, null, persistent);
		 sectionPoint.addSectionRowCol(sectionRowF, sectionColF, this);
//rowColFix
//					if(rowIsNext)
//						sectionPoint.setRowOrColIndexF(null, ROW, gridRowF);
//					if(colIsNext)
//						sectionPoint.setRowOrColIndexF(null, COL, gridColF);
//
		 return sectionPoint;
		}
	   }
	   return null;
	  }
	  catch(Exception e)
	  {
	   StsException.outputException("StsSection.getNextGridCrossing() failed.",
		e, StsException.WARNING);
	   return null;
	  }
	 }
	 */

	/** Add sectionRowCrossings and XY gridCrossings */
	/*
	 private void addRowSectionPoints(StsGridSectionPoint sectionPoint0, StsGridSectionPoint sectionPoint1,
			  StsEdgeLinkable edge, StsXYSurfaceGridable grid, StsObjectList sectionPoints,
			  int extrapolate, boolean persistent)
	 {
	  StsGridSectionPoint sectionPoint;
	  StsList gridSectionPoints;
	  int rowStart, rowEnd;
	  int dRow;
	  int sectionRowConnect = NONE;

	  float rowF0 = sectionPoint0.getRowF(this);
	  float rowF1 = sectionPoint1.getRowF(this);

	  if(rowF0 < 0 || rowF0 > nRows-1 || rowF1 < 0 || rowF1 > nRows-1) return;

	  if(rowF1 == rowF0)
	  {
	   gridSectionPoints = StsGridSectionPoint.getGridCrossings(grid, this,
		sectionPoint0, sectionPoint1, edge, extrapolate, persistent, MINUS);
	   addSectionPointsToList(gridSectionPoints, sectionPoints);
	   return;
	  }
	  else if(rowF1 > rowF0)
	  {
	   rowStart = StsMath.ceiling(rowF0);
	   rowEnd = StsMath.floor(rowF1);
	   dRow = 1;
	   if(rowStart > rowEnd)
	   {
		gridSectionPoints = StsGridSectionPoint.getGridCrossings(grid, this,
		 sectionPoint0, sectionPoint1, edge, extrapolate, persistent, MINUS);
		addSectionPointsToList(gridSectionPoints, sectionPoints);
		return;
	   }
	   sectionRowConnect = MINUS;
	  }
	  else
	  {
	   rowStart = StsMath.floor(rowF0);
	   rowEnd = StsMath.ceiling(rowF1);
	   dRow = -1;
	   if(rowStart < rowEnd)
	   {
		gridSectionPoints = StsGridSectionPoint.getGridCrossings(grid, this,
		 sectionPoint0, sectionPoint1, edge, extrapolate, persistent, MINUS);
		addSectionPointsToList(gridSectionPoints, sectionPoints);
		return;
	   }
	   sectionRowConnect = PLUS;
	  }

	  float colf0 = sectionPoint0.getColF(this);
	  float colf1 = sectionPoint1.getColF(this);
	  StsGridSectionPoint lastSectionPoint = sectionPoint0;

	  int row = rowStart;
	  while(true)
	  {
	   float ff = (row - rowF0)/(rowF1 - rowF0);
	   if(ff > 1.0f) break;
	 sectionPoint = StsGridSectionPoint.sectionInterpolate(sectionPoint0, sectionPoint1, (double)ff, edge, this, persistent);
	   sectionPoint.setRowOrColIndex(this, ROW, row, false);
	   gridSectionPoints = StsGridSectionPoint.getGridCrossings(grid, this, lastSectionPoint, sectionPoint, edge, extrapolate, persistent, MINUS);
	   addSectionPointsToList(gridSectionPoints, sectionPoints);
	   sectionPoints.add(sectionPoint);
	   lastSectionPoint = sectionPoint;
	   row += dRow;
	  }

//        sectionPoint1.initializeGridEdgePoint();
//        sectionPoint1.setGridRowColF();
	  gridSectionPoints = StsGridSectionPoint.getGridCrossings(grid, this,
	   lastSectionPoint, sectionPoint1, edge, extrapolate, persistent, MINUS);
	  addSectionPointsToList(gridSectionPoints, sectionPoints);
	 }
	 */
	/*
	 private void checkAddSectionPoint(StsObjectList sectionPoints, StsGridSectionPoint newPoint)
	 {
	  StsGridSectionPoint lastPoint = (StsGridSectionPoint)sectionPoints.getLast();
	  if(newPoint != lastPoint) sectionPoints.add(newPoint);
	 }
	 */
	private void getGridCrossings(StsXYSurfaceGridable grid,
								  StsGridSectionPoint sectionPoint0, StsGridSectionPoint sectionPoint1,
								  StsEdgeLinkable edge, int extrapolate, boolean persistent, int endFlags,
								  StsObjectList sectionPoints)
	{
		StsList gridSectionPoints = StsGridSectionPoint.getGridCrossings(grid, this,
			sectionPoint0, sectionPoint1, edge, extrapolate, persistent, endFlags);
		addSectionPointsToList(gridSectionPoints, sectionPoints);
	}

	private void checkRemoveDyingFaultPoints(StsSectionEdge edge, StsObjectList sectionPoints)
	{
		StsSurfaceVertex vertex;
		StsGridSectionPoint edgePoint;

		vertex = edge.getPrevVertex();
		if (vertex.onDyingFault())
		{
			edgePoint = (StsGridSectionPoint) sectionPoints.getFirst();
			while (edgePoint != null)
			{
				if (edgePoint.isRowOrCol(null))
				{
					break;
				}
				sectionPoints.delete(0);
				edgePoint = (StsGridSectionPoint) sectionPoints.getFirst();
			}
		}

		vertex = edge.getNextVertex();
		if (vertex.onDyingFault())
		{
			edgePoint = (StsGridSectionPoint) sectionPoints.getLast();
			while (edgePoint != null)
			{
				if (edgePoint.isRowOrCol(null))
				{
					break;
				}
				sectionPoints.deleteLast();
				edgePoint = (StsGridSectionPoint) sectionPoints.getLast();
			}
		}
	}

	private StsObjectList constructRowColSurfaceEdgePoints(StsEdge edge, StsXYSurfaceGridable grid, boolean persistent)
	{
		StsGridSectionPoint surfacePoint0, surfacePoint1;
		StsGridPoint gridPoint;
		StsGridSectionPoint edgePoint;
		StsObjectList edgePoints = null;
		int nPnts;
		int connect;
		float i0F = -1, i1F = -1, j0F = -1, j1F = -1;
		int i, j, i0, i1, j0, j1;

		try
		{
			int rowOrCol = geometry.getRowOrCol();

			if (rowOrCol == ROW)
			{
				surfacePoint0 = edge.getPrevVertex().getSurfacePoint();
				i0F = surfacePoint0.getRowF(null);
				j0F = surfacePoint0.getColF(null);
//rowColFix
				/*
					i0F = StsMath.roundOffInteger(i0F);
					surfacePoint0.setRowOrColIndexF(null, ROW, i0F);
					if(StsMath.isIntegral(j0F))
					{
					 j0F = StsMath.roundOffInteger(j0F);
					 surfacePoint0.setRowOrColIndexF(null, COL, j0F);
					}
				 */
				surfacePoint1 = edge.getNextVertex().getSurfacePoint();
				i1F = surfacePoint1.getRowF(null);
//rowColFix
				/*
					i1F = StsMath.roundOffInteger(i1F);
					surfacePoint1.setRowOrColIndexF(null, ROW, i1F);
					if(i1F != i0F)
					{
					 StsException.systemError("StsSection.constructRowColSurfaceEdgePoints() failed." +
				 "SurfaceEdge: " + edge.getLabel() + " is not on a row, but is from row: " + i0F + " to: " + i1F);
					}
				 */
				i = (int) i0F;

				j1F = surfacePoint1.getColF(null);
//rowColFix
				/*
					if(StsMath.isIntegral(j1F))
					{
					 j1F = StsMath.roundOffInteger(j1F);
					 surfacePoint1.setRowOrColIndexF(null, COL, j1F);
					}
				 */
				if (j1F >= j0F)
				{
					j0 = StsMath.ceiling(j0F);
					j1 = StsMath.floor(j1F);

					//                connect = MINUS;
					nPnts = j1 - j0 + 1;

					/** Allocate space for grid points plus possibly two end points. */
					edgePoints = new StsObjectList(nPnts + 2, 2);

					/** If point is on row and col, set direction orthogonal to edge;
					 *  otherwise leave its direction same as this edge */
					if (surfacePoint0.isRowAndCol(null))
					{
						j0++;
					}
					if (surfacePoint1.isRowAndCol(null))
					{
						j1--;

					}
					edgePoints.add(surfacePoint0);

					for (j = j0; j <= j1; j++)
					{
						edgePoint = getGapOrGridPoint(i, j, edge, grid, persistent);
						edgePoints.add(edgePoint);
					}

					edgePoints.add(surfacePoint1);
				}
				else
				{
					j0 = StsMath.floor(j0F);
					j1 = StsMath.ceiling(j1F);

					nPnts = j0 - j1 + 1;

					edgePoints = new StsObjectList(nPnts + 2, 2);

					if (surfacePoint0.isRowAndCol(null))
					{
						j0--;
					}
					if (surfacePoint1.isRowAndCol(null))
					{
						j1++;

					}
					edgePoints.add(surfacePoint0);

					for (j = j0; j >= j1; j--)
					{
						edgePoint = getGapOrGridPoint(i, j, edge, grid, persistent);
						edgePoints.add(edgePoint);
					}
					edgePoints.add(surfacePoint1);
				}
			}
			else if (rowOrCol == COL)
			{
				surfacePoint0 = edge.getPrevVertex().getSurfacePoint();
				j0F = surfacePoint0.getColF(null);
				i0F = surfacePoint0.getRowF(null);
//rowColFix
				/*
					j0F = StsMath.roundOffInteger(j0F);
					surfacePoint0.setRowOrColIndexF(null, COL, j0F);
					if(StsMath.isIntegral(i0F))
					{
					 i0F = StsMath.roundOffInteger(i0F);
					 surfacePoint0.setRowOrColIndexF(null, ROW, i0F);
					}
				 */
				surfacePoint1 = edge.getNextVertex().getSurfacePoint();
				j1F = surfacePoint1.getColF(null);
				i1F = surfacePoint1.getRowF(null);
//rowColFix
				/*
					j1F = StsMath.roundOffInteger(j1F);
					surfacePoint1.setRowOrColIndexF(null, COL, j1F);
					if(StsMath.isIntegral(i1F))
					{
					 i1F = StsMath.roundOffInteger(i1F);
					 surfacePoint1.setRowOrColIndexF(null, ROW, i1F);
					}
				 */
				if (j1F != j0F)
				{
					StsException.systemError("StsSection.constructRowColSurfaceEdgePoints() failed." +
											 "SurfaceEdge: " + edge.getLabel() + " is not on a col, but is from col: " +
											 j0F + " to: " + j1F);
				}

				j = (int) j0F;

				if (i1F >= i0F)
				{
					i0 = StsMath.ceiling(i0F);
					i1 = StsMath.floor(i1F);
					nPnts = i1 - i0 + 1;

					edgePoints = new StsObjectList(nPnts + 2, 2);

					if (surfacePoint0.isRowAndCol(null))
					{
						i0++;
					}
					if (surfacePoint1.isRowAndCol(null))
					{
						i1--;

					}
					edgePoints.add(surfacePoint0);

					for (i = i0; i <= i1; i++)
					{
						edgePoint = getGapOrGridPoint(i, j, edge, grid, persistent);
						edgePoints.add(edgePoint);
					}

					edgePoints.add(surfacePoint1);
				}
				else
				{
					i0 = StsMath.floor(i0F);
					i1 = StsMath.ceiling(i1F);
					nPnts = i0 - i1 + 1;

					edgePoints = new StsObjectList(nPnts + 2, 2);

					if (surfacePoint0.isRowAndCol(null))
					{
						i0--;
					}
					if (surfacePoint1.isRowAndCol(null))
					{
						i1++;

					}
					edgePoints.add(surfacePoint0);

					for (i = i0; i >= i1; i--)
					{
						edgePoint = getGapOrGridPoint(i, j, edge, grid, persistent);
						edgePoints.add(edgePoint);
					}

					edgePoints.add(surfacePoint1);
				}
			}
			else
			{
				StsException.systemError("StsSection.constructRowColSurfaceEdgePoints() failed.\n" +
										 " Edge is not aligned with row or col:\n" +
										 " SurfaceEdge: " + edge.getLabel() + i0F + "," + j0F + " to: " + i1F + "," +
										 j1F);
			}

			edgePoints.trimToSize();
			return edgePoints;
		}
		catch (Exception e)
		{
			StsException.outputException("StsSection.constructRowColSurfaceEdgePoints() failed.\n" +
										 " Edge is not aligned with row or col:\n" +
										 " SurfaceEdge: " + edge.getLabel() + i0F + "," + j0F + " to: " + i1F + "," +
										 j1F,
										 e, StsException.WARNING);
			return null;
		}
	}

	/** We may have 1 or 2 grids; try to get a gap value from these; failing
	 *  that get it from the surface grid.
	 */
	private StsGridSectionPoint getGapOrGridPoint(int row, int col, StsEdge edge, StsXYSurfaceGridable grid, boolean persistent)
	{
		StsGridSectionPoint edgePoint;
		StsPoint point;

		try
		{
			if (grid instanceof StsBlockGrid)
			{
				StsBlockGrid blockGrid = (StsBlockGrid) grid;
				point = blockGrid.getComputePoint(row, col, true);
				if (point != null)
				{
					edgePoint = new StsGridSectionPoint(point, row, col, edge, blockGrid, persistent);
					edgePoint.addSectionRowCol(this);
					return edgePoint;
				}
			}
			// grid instanceof StsGrid || point == null (from above)
//			StsModelSurface surface = (StsModelSurface) grid;
            point = grid.getPoint(row, col);
			edgePoint = new StsGridSectionPoint(point, row, col, edge, grid, persistent);
			edgePoint.addSectionRowCol(this);
			return edgePoint;
		}
		catch (Exception e)
		{
			StsException.outputException("StsSection.getGapOrGridPoint() failed",
										 e, StsException.WARNING);
			return null;
		}
	}

	/*
	 private static StsObjectRefList convertListToRefList(StsObjectList edgePointsList)
	 {
	  int n;

	  int nPnts = edgePointsList.getSize();
	  for(n = 0; n < nPnts; n++)
	  {
	   StsObject edgePoint = (StsObject)edgePointsList.getElement(n);
	   if(!edgePoint.isPersistent()) currentModel.add(edgePoint);
	  }

	  StsObjectRefList edgePoints = StsObjectRefList.constructor(edgePointsList);
	  return edgePoints;
	 }
	 */
	private void addSectionPointsToList(StsList gridSectionPoints, StsObjectList sectionPoints)
	{
		StsObject point;

		if (gridSectionPoints == null || sectionPoints == null)
		{
			return;
		}

		StsObject lastPoint = (StsObject) sectionPoints.getLast();
		int n = 0;

		if (lastPoint != null)
		{
			while ( (point = (StsObject) gridSectionPoints.getElement(n)) == lastPoint)
			{
				n++;
			}
		}
		int nPnts = gridSectionPoints.getSize();
		for (; n < nPnts; n++)
		{
			point = (StsObject) gridSectionPoints.getElement(n);
			sectionPoints.add(point);
		}
	}

	private void addSectionPointToList(StsGridSectionPoint gridSectionPoint, StsObjectList sectionPoints)
	{
		if (gridSectionPoint == null || sectionPoints == null)
		{
			return;
		}
		StsObject lastPoint = (StsObject) sectionPoints.getLast();
		if (lastPoint != null && gridSectionPoint == lastPoint)
		{
			return;
		}
		sectionPoints.add(gridSectionPoint);
	}

// Methods for handling boundary sections

	static public boolean clipEdgeToBoundaries(StsSectionEdge edge)
	{
		boolean ok;
		StsSectionEdge newEdge;
		boolean deleteEdge = false;

		if (currentModel == null)
		{
			return false;
		}
		StsClass sections = currentModel.getCreateStsClass(StsSection.class);
		int nSections = sections.getSize();

		StsList intersections = new StsList(4, 2);

		for (int i = 0; i < nSections; i++)
		{
			StsSection section = (StsSection) sections.getElement(i);
			if (section.getType() == BOUNDARY)
			{
				ok = section.clipEdgeToBoundary(edge, intersections);
				if (!ok)
				{
					return false;
				}
			}
		}

		StsList edgePoints = edge.getEdgePointsList();
		int nEdgePoints = edgePoints.getSize();

		int nIntersections = intersections.getSize();

		StsSurfaceVertex startVertex = null;
		int startIndex = 0;
		for (int n = 0; n < nIntersections; n++)
		{
			Intersection intersection = (Intersection) intersections.getElement(n);
			if (intersection.direction == OUTSIDE)
			{
				newEdge = edge.makeSubEdge(startVertex, startIndex, intersection.vertex, intersection.vertexIndex);
				if (newEdge != null)
				{
					deleteEdge = true;
				}
				startVertex = null;
			}
			else if (intersection.direction == INSIDE)
			{
				startVertex = intersection.vertex;
				startIndex = intersection.vertexIndex + 1;
			}
		}

		if (startVertex != null)
		{
			newEdge = edge.makeSubEdge(startVertex, startIndex, null, nEdgePoints - 1);
			if (newEdge != null)
			{
				deleteEdge = true;
			}
		}

		if (deleteEdge)
		{
			edge.deleteEdgeNotPoints();
//        else
//            edge.constructPoints();

		}
		return true;
	}

	public boolean clipEdgeToBoundary(StsSectionEdge edge, StsList intersections)
	{
		StsObjectList edgePoints;
		StsSurfaceVertex firstVertex, lastVertex;
		int lastSide, side;
		StsSectionPoint lastSectionPoint, sectionPoint;
		Intersection intersection;
		boolean intersected;

		try
		{
			StsPoint[] points = edge.getPoints();
			if (points == null)
			{
				return false;
			}
			int nPoints = points.length;

			sectionPoint = new StsSectionPoint(points[0]);
			computeNearestPoint(sectionPoint);
			side = getSide(sectionPoint);

			for (int i = 1; i < nPoints; i++)
			{
				lastSectionPoint = sectionPoint;
				lastSide = side;
				sectionPoint = new StsSectionPoint(points[i]);
				sectionPoint.sectionColF = lastSectionPoint.sectionColF;
				computeNearestPoint(sectionPoint);

				int direction = getDirection(i, nPoints, lastSectionPoint, sectionPoint);

				if (direction != NONE)
				{
					if (i == 1)
					{
						firstVertex = edge.getPrevVertex();
						if (firstVertex != null && firstVertex.isConnectedToSection(this))
						{
							continue;
						}
					}
					if (i == nPoints - 1)
					{
						lastVertex = edge.getNextVertex();
						if (lastVertex != null && lastVertex.isConnectedToSection(this))
						{
							continue;
						}
					}
					addIntersection(edge, lastSectionPoint, sectionPoint, i, direction, intersections);
				}
			}

			return true;
		}
		catch (Exception e)
		{
			StsException.outputException(e, StsException.WARNING);
			return false;
		}
	}

	private int getSide(StsSectionPoint sectionPoint)
	{
		if (sectionPoint.side == RIGHT)
		{
			return INSIDE;
		}
		else
		{
			return OUTSIDE;
		}
	}

	private int getDirection(int index, int nPnts, StsSectionPoint lastSectionPoint, StsSectionPoint sectionPoint)
	{
		int side = sectionPoint.side;
		int lastSide = lastSectionPoint.side;

		if (side != lastSide)
		{
			if (side == RIGHT)
			{
				return INSIDE;
			}
			else
			{
				return OUTSIDE;
			}
		}

		float avgGridInc = currentModel.getGridDefinition().getAvgInc();
		float distance = sectionPoint.distance;
		float lastDistance = lastSectionPoint.distance;

		if (index == 1 && lastDistance < avgGridInc && lastDistance <= distance)
		{
			if (lastSide == RIGHT)
			{
				return INSIDE;
			}
			else
			{
				return OUTSIDE;
			}
		}
		else if (index == nPnts - 1 && distance < avgGridInc && distance <= lastDistance)
		{
			if (side == RIGHT)
			{
				return OUTSIDE;
			}
			else
			{
				return INSIDE;
			}
		}
		else
		{
			return NONE;
		}
	}

	private boolean addEndIntersection(StsSectionEdge edge, StsSectionPoint sectionPoint,
									   StsSurfaceVertex vertex, int index, int direction, StsList intersections)
	{
		Intersection intersection;

		if (sectionPoint.offSection)
		{
			return false;
		}
		intersection = new Intersection(edge, (float) index, direction, vertex);
		intersection.addToList(intersections);
		return true;
	}

	public void setIsVisible(boolean b)
	{
		isVisible = b;
        checkSetSectionClass();
        if (b && currentModel != null)
		{
			switch (type)
			{
				case BOUNDARY:
					sectionClass.setDisplayBoundarySections(true);
					break;
				case FAULT:
					sectionClass.setDisplayFaultSections(true);
					break;
				case FRACTURE:
					sectionClass.setDisplayFractureSections(true);
					break;					
				case AUXILIARY:
					sectionClass.setDisplayAuxiliarySections(true);
					break;
			}
		}
		toggleLines(b);
	}

	public boolean getIsVisible()
	{
		return isVisible;
	}

	private boolean addIntersection(StsSectionEdge edge, StsSectionPoint sectionPoint0, StsSectionPoint sectionPoint1,
									int index1, int direction, StsList intersections)
	{
		StsSurfaceVertex vertex0, vertex1, vertex;

		Intersection intersection;

		float d0 = sectionPoint0.distance * sectionPoint0.side;
		float d1 = sectionPoint1.distance * sectionPoint1.side;
		if (d0 != 0.0f && d0 == d1)
		{
			return false; /** Intersection failed: line is parallel to section */
		}

		if (sectionPoint0.offSection && sectionPoint1.offSection)
		{
			return false;
		}

		StsSectionPoint sectionPoint = new StsSectionPoint();
//        sectionPoint.point = new StsPoint();

		float f, lastf;

		if (Math.abs(d0) < StsParameters.roundOff)
		{
			f = 0.0f;
		}
		else if (Math.abs(d1) < StsParameters.roundOff)
		{
			f = 1.0f;
		}
		else
		{
			f = -d0 / (d1 - d0);

		}

		float df = 1.0f;
		while (df > 0.001f)
		{
			sectionPoint.point.interpolatePoints(sectionPoint0.point, sectionPoint1.point, f);
			computeNearestPoint(sectionPoint);
			lastf = f;
			f = sectionPoint.nearestPoint.projectedInterpolationFactor(3, sectionPoint0.point, sectionPoint1.point);
			df = Math.abs(f - lastf);
		}

		if (sectionPoint.offSection)
		{
			return false;
		}

		f = StsMath.minMax(f, 0.001f, 0.999f);

		float sectionColF = sectionPoint.sectionColF;
		StsGridSectionPoint gridPoint = new StsGridSectionPoint(sectionPoint.point, this, true);
		gridPoint.setRowOrColIndexF(this, COL, sectionColF);
		vertex = new StsSurfaceVertex(gridPoint);

		int index0 = index1 - 1;
		intersection = new Intersection(edge, index0 + f, direction, vertex);
		intersection.addToList(intersections);
		return true;
	}

	/*
	 private boolean addIntersection(StsSectionEdge edge, StsSectionPoint sectionPoint0, StsSectionPoint sectionPoint1,
			 int index1, int direction, StsList intersections)
	 {
	  StsSurfaceVertex vertex0, vertex1, vertex;

	  Intersection intersection;

	  float d0 = sectionPoint0.distance*sectionPoint0.side;
	  float d1 = sectionPoint1.distance*sectionPoint1.side;
	  if(d0 != 0.0f && d0 == d1) return false; // Intersection failed: line is parallel to section

	  int index0 = index1-1;
	  vertex0 = edge.getVertexAtPoint(index0);
	  vertex1 = edge.getVertexAtPoint(index1);

	  if(vertex0 != null && vertex0.getSectionEdge() == sectionEdge)
	  {
	   if(sectionPoint0.offSection) return false;
	   intersection = new Intersection(edge, (float)index0, direction, vertex0);
	   intersection.addToList(intersections);
	   return true;
	  }
	  else if(vertex1 != null && vertex1.getSectionEdge() == sectionEdge)
	  {
	   if(sectionPoint1.offSection) return false;
	   intersection = new Intersection(edge, (float)(index1), direction, vertex1);
	   intersection.addToList(intersections);
	   return true;
	  }

	  if(sectionPoint0.offSection && sectionPoint1.offSection) return false;

	  StsSectionPoint sectionPoint = new StsSectionPoint();
	  sectionPoint.point = new StsPoint();

	  float f, lastf;

	  if(Math.abs(d0) < StsParameters.roundOff)
	   f = 0.0f;
	  else if(Math.abs(d1) < StsParameters.roundOff)
	   f = 1.0f;
	  else
	   f = -d0/(d1 - d0);


	  float df = 1.0f;
	  while(df > 0.001f)
	  {
	   sectionPoint.point.interpolatePoints(sectionPoint0.point, sectionPoint1.point, f);
	   computeNearestPoint(sectionPoint);
	   lastf = f;
	 f = sectionPoint.nearestPoint.projectedInterpolationFactor(3, sectionPoint0.point, sectionPoint1.point);
	   df = Math.abs(f - lastf);
	  }

	  if(sectionPoint.offSection) return false;

	  if(vertex0 != null && f <= 0.0f && vertex0.isFullyConnected()) return false;
	  if(vertex1 != null && f >= 1.0f && vertex1.isFullyConnected()) return false;

	  float sectionColF = sectionPoint.sectionColF;
	  vertex = new StsSurfaceVertex(sectionPoint.point, sectionEdge, sectionColF);

	  intersection = new Intersection(edge, index0+f, direction, vertex);
	  intersection.addToList(intersections);
	  return true;
	 }
	 */
	/*
	 private boolean addIntersection(StsSectionEdge edge, StsSectionPoint sectionPoint0, StsSectionPoint sectionPoint1,
			 int index1, int direction, StsList intersections)
	 {
	  StsSurfaceVertex vertex0, vertex1, vertex;

	  Intersection intersection;

	  float d0 = sectionPoint0.distance*sectionPoint0.side;
	  float d1 = sectionPoint1.distance*sectionPoint1.side;
	  if(d0 != 0.0f && d0 == d1) return false; // Intersection failed: line is parallel to section

	  int index0 = index1-1;
	  vertex0 = edge.getVertexAtPoint(index0);
	  if(vertex0 != null && vertex0.getSectionEdge() == sectionEdge) return false; // already connected

	  vertex1 = edge.getVertexAtPoint(index1);
	  if(vertex1 != null && vertex1.getSectionEdge() == sectionEdge) return false; // already connected

	  if(vertex0 != null)
	  {
	   if(sectionPoint0.offSection) return false;
	   intersection = new Intersection(edge, (float)index0, direction, vertex0);
	   intersection.addToList(intersections);
	   return true;
	  }

	  if(vertex1 != null)
	  {
	   if(sectionPoint1.offSection) return false;
	   intersection = new Intersection(edge, (float)(index1), direction, vertex1);
	   intersection.addToList(intersections);
	   return true;
	  }

	  if(sectionPoint0.offSection && sectionPoint1.offSection) return false;

	  StsSectionPoint sectionPoint = new StsSectionPoint();
	  sectionPoint.point = new StsPoint();

	  float f, lastf;

	  if(Math.abs(d0) < StsParameters.roundOff)
	   f = 0.0f;
	  else if(Math.abs(d1) < StsParameters.roundOff)
	   f = 1.0f;
	  else
	   f = -d0/(d1 - d0);


	  float df = 1.0f;
	  while(df > 0.001f)
	  {
	   sectionPoint.point.interpolatePoints(sectionPoint0.point, sectionPoint1.point, f);
	   computeNearestPoint(sectionPoint);
	   lastf = f;
	 f = sectionPoint.nearestPoint.projectedInterpolationFactor(3, sectionPoint0.point, sectionPoint1.point);
	   df = Math.abs(f - lastf);
	  }

	  if(sectionPoint.offSection) return false;

	  if(vertex0 != null && f <= 0.0f && vertex0.isFullyConnected()) return false;
	  if(vertex1 != null && f >= 1.0f && vertex1.isFullyConnected()) return false;

	  float sectionColF = sectionPoint.sectionColF;
	  vertex = new StsSurfaceVertex(sectionPoint.point, sectionEdge, sectionColF);

	  intersection = new Intersection(edge, index0+f, direction, vertex);
	  intersection.addToList(intersections);
	  return true;
	 }
	 */
	class Intersection
	{
		float sectionColF;
		int direction;
		StsSurfaceVertex vertex;
		int vertexIndex;

		Intersection(StsSectionEdge intersectingEdge, float sectionColF, int direction, StsSurfaceVertex vertex)
		{
			this.sectionColF = sectionColF;
			this.direction = direction;
			this.vertex = vertex;
			vertexIndex = (int) sectionColF / intersectingEdge.getNPointsPerVertex();
			StsXYSurfaceGridable surface = intersectingEdge.getSurface();
			if (surface != null)
			{
				StsSectionEdge intersectedEdge = getSurfaceSectionEdge(surface);
				vertex.addEdgeAssociation(intersectedEdge);
			}
		}

		void addToList(StsList intersections)
		{
			int nIntersections = intersections.getSize();
			for (int n = 0; n < nIntersections; n++)
			{
				Intersection intersection = (Intersection) intersections.getElement(n);
				if (intersection.sectionColF > sectionColF)
				{
					intersections.insertBefore(n, this);
					return;
				}
			}
			intersections.add(this);
		}
	}

	public void constructBlockSides()
	{
		constructBlockSides(RIGHT);
		if (!isBoundary())
		{
			constructBlockSides(LEFT);
		}
	}

	private void constructBlockSides(int side)
	{
		StsLine line, nextLine;
		StsBlockSide blockSide;
		int n;

		StsObjectList sectionLines = getOrderedSectionLines(side);
		if (sectionLines == null)
		{
			return;
		}
		int nSectionLines = sectionLines.getSize();

		StsObjectRefList blockSides;

		if (side == LEFT)
		{
			if (leftBlockSides != null)
			{
				leftBlockSides.deleteAll();
			}
			leftBlockSides = StsObjectRefList.constructor(nSectionLines - 1, 1, "leftBlockSides", this);
			blockSides = leftBlockSides;
		}
		else
		{
			if (rightBlockSides != null)
			{
				rightBlockSides.deleteAll();
			}
			rightBlockSides = StsObjectRefList.constructor(nSectionLines - 1, 1, "rightBlockSides", this);
			blockSides = rightBlockSides;
		}

		nextLine = (StsLine) sectionLines.getElement(0);
		for (n = 1; n < nSectionLines; n++)
		{
			line = nextLine;
			nextLine = (StsLine) sectionLines.getElement(n);
			blockSide = new StsBlockSide(line, nextLine, this, side);
			blockSides.add(blockSide);
		}
	}

	private StsObjectList getOrderedSectionLines(int side)
	{
		StsLine line, orderedLine;

		StsObjectRefList sectionLines;
		if (side == RIGHT)
		{
			sectionLines = rightLines;
		}
		else if (side == LEFT)
		{
			sectionLines = leftLines;
		}
		else
		{
			return null;
		}

		int nSectionLines = sectionLines.getSize();
		StsObjectList orderedLines = new StsObjectList(nSectionLines + 2);

		if (side == RIGHT)
		{
			orderedLines.add(getFirstLine());
		}
		else
		{
			orderedLines.add(getLastLine());

		}
		for (int n = 0; n < nSectionLines; n++)
		{
			line = (StsLine) sectionLines.getElement(n);
			orderedLines.add(line);
		}

		if (side == RIGHT)
		{
			orderedLines.add(getLastLine());
		}
		else
		{
			orderedLines.add(getFirstLine());

		}
		return orderedLines;
	}

	public StsBlockSide getNextBlockSide(StsLine line)
	{
		StsObjectRefList blockSides;

		if (line == getFirstLine())
		{
			return (StsBlockSide) rightBlockSides.getFirst();
		}
		else if (line == getLastLine())
		{
			return (StsBlockSide) leftBlockSides.getFirst();
		}
		else
		{
			if (line.getOnSection() != this)
			{
				return null;
			}
			int side = line.getSectionSide();
			if (side == RIGHT)
			{
				blockSides = rightBlockSides;
			}
			else if (side == LEFT)
			{
				blockSides = leftBlockSides;
			}
			else
			{
				return null;
			}

			int nBlockSides = blockSides.getSize();

			for (int n = 0; n < nBlockSides; n++)
			{
				StsBlockSide blockSide = (StsBlockSide) blockSides.getElement(n);
				if (blockSide.getPrevLine() == line)
				{
					return blockSide;
				}
			}
			StsException.systemError("StsSection.getNextBlockSide() failed." +
									 " Couldn't find next blockSide from blockSide: " + line.getLabel());
			return null;
		}
	}

	public StsBlockSide getPrevBlockSide(StsLine line)
	{
		StsObjectRefList blockSides;

		if (line == getFirstLine())
		{
			return (StsBlockSide) leftBlockSides.getLast();
		}
		else if (line == getLastLine())
		{
			return (StsBlockSide) rightBlockSides.getLast();
		}
		else
		{
			if (line.getOnSection() != this)
			{
				return null;
			}
			int side = line.getSectionSide();
			if (side == RIGHT)
			{
				blockSides = rightBlockSides;
			}
			else if (side == LEFT)
			{
				blockSides = leftBlockSides;
			}
			else
			{
				return null;
			}

			int nBlockSides = blockSides.getSize();

			for (int n = 0; n < nBlockSides; n++)
			{
				StsBlockSide blockSide = (StsBlockSide) blockSides.getElement(n);
				if (blockSide.getNextLine() == line)
				{
					return blockSide;
				}
			}
			StsException.systemError("StsSection.getPrevBlockSide() failed." +
									 " Couldn't find next blockSide from blockSide: " + line.getLabel());
			return null;
		}
	}

	public void addZoneSide(StsZoneSide zoneSide)
	{
		if (zoneSide.getSide() == RIGHT)
		{
			if (rightZoneSides == null)
			{
				rightZoneSides = new StsList(5, 5);
			}
			rightZoneSides.add(zoneSide);
		}
		else
		{
			if (leftZoneSides == null)
			{
				leftZoneSides = new StsList(5, 5);
			}
			leftZoneSides.add(zoneSide);
		}
	}
/*
    public void addBlockGridPolygons()
    {
        byte geometryType = getGeometryType();
        addBlockGridPolygons(leftZoneSides, geometryType);
        addBlockGridPolygons(rightZoneSides, geometryType);
    }

    private void addBlockGridPolygons(StsObjectRefList zoneSides, byte geometryType)
    {
        for(int n = 0; n < zoneSides.getSize(); n++)
        {
            StsZoneSide zoneSide = (StsZoneSide)zoneSides.getElement(n);
            zoneSide.addBlockGridPolygons(geometryType);
        }
    }
*/
    public int getEnd(StsLine line)
	{
		if (line == getFirstLine())
		{
			return MINUS;
		}
		else if (line == getLastLine())
		{
			return PLUS;
		}
		else
		{
			return NONE;
		}
	}

	private void setEdgeEndCols(StsEdge edge)
	{
		if (edge == null)return;
		edge.getFirstEdgePoint().addSectionRowCols();
		edge.getLastEdgePoint().addSectionRowCols();
	}

    /*
	 private void constructIntersectedEdges()
	 {
	  StsSurface[] surfaces = model.getSurfaces(StsModelSurface.MODEL, true);
	  int nSurfaces = surfaces.length;
	  intersectedEdges = new StsSectionEdge[nSurfaces];

	  for(int n = 0; n < nSurfaces; n++)
	   intersectedEdges[n] = new StsSectionEdge(this, surfaces[n], false);
	 }
	 */
// 3D display routines

	public void display(StsGLPanel3d glPanel3d, boolean displaySectionEdges)
	{
		if (!isVisible || !displayOK)
        {
            checkDeleteDisplays();
            return;
        }
		byte projectZDomain = currentModel.getProject().getZDomain();
		if (projectZDomain != zDomainDisplayed)
		{
			checkDeleteDisplays();
            if(!canDraw()) return;
            zDomainDisplayed = projectZDomain;
            constructSection();
		}
		useDisplayLists = currentModel.useDisplayLists;
		if (!useDisplayLists && usingDisplayLists)
		{
			deleteDisplayLists();
			usingDisplayLists = false;
		}

		GL gl = glPanel3d.getGL();
		boolean isDepth = projectZDomain == StsProject.TD_DEPTH;
		if (sectionPatch != null && drawSurface)
		{
            if(StsGLPanel.debugProjectionMatrix) StsException.systemDebug("viewShift called for displaySurface of " + toString());
            glPanel3d.setViewShift(gl, StsGraphicParameters.sectionShift);
			sectionPatch.displaySurface(currentModel, displaySide, isDepth);
			glPanel3d.resetViewShift(gl);
		}

		if (sectionEdges != null && (displaySectionEdges || drawEdges))
		{
			int nSectionEdges = sectionEdges.getSize();
			for (int n = 0; n < nSectionEdges; n++)
			{
				StsSectionEdge sectionEdge = (StsSectionEdge) sectionEdges.getElement(n);
				StsXYSurfaceGridable surface = sectionEdge.getSurface();
				if(surface != null && surface.getIsVisible())
				{
					sectionEdge.display(glPanel3d);
				}
			}
		}
		if (sectionPatch != null && drawGrid)
		{
			sectionPatch.displayGrid(glPanel3d);
		}
		boolean debug = currentModel.getBooleanProperty("debugDisplayZoneSides");
		if (debug)
		{
			StsSurfaceEdge currentEdge = StsSurfaceEdge.getCurrentEdge();
			if (currentEdge != null)
			{
				StsZoneSide zoneSide = currentEdge.getAZoneSide();
				if (zoneSide != null && zoneSide.getSide() == displaySide)
				{
					zoneSide.debugDisplay(currentModel.win3d);
				}
			}
		}
	}

    protected void checkDeleteDisplays()
	{
		//       if(surfaceChanged || textureChanged) return;
		deleteDisplayLists();
		//       deleteTextureTileSurface(gl);
	}
/*
    private void checkAdjustTimeOrDepth()
    {
        StsSeismicVelocityModel velocityModel = currentModel.project.getSeismicVelocityModel();
        if(velocityModel != null) adjustTimeOrDepth(velocityModel);
    }
*/
    public boolean adjustTimeOrDepth(StsSeismicVelocityModel velocityModel)
	{
        boolean isOriginalDepth = isOriginalDepth();
        StsLine line = getFirstLine();
		line.adjustTimeOrDepth(velocityModel, isOriginalDepth);
		line = getLastLine();
		line.adjustTimeOrDepth(velocityModel, isOriginalDepth);
		int nSectionEdges = sectionEdges.getSize();
		for (int n = 0; n < nSectionEdges; n++)
		{
			StsSectionEdge sectionEdge = (StsSectionEdge) sectionEdges.getElement(n);
			sectionEdge.adjustTimeOrDepth(velocityModel, isOriginalDepth);
		}
		return true;
	}

    public byte getGeometryType()
    {
        return geometry.getGeometryType();
    }

    public StsGridBoundingBox getGridBoundingBox()
    {
        return new StsGridBoundingBox(0, nRows-1, 0, nCols-1);
    }

    public StsFieldBean[] getDisplayFields()
	{
		if (type != StsSection.FAULT)
		{
            if(displayFields != null) return displayFields;
            displayFields = new StsFieldBean[]
		    {
                new StsBooleanFieldBean(StsSection.class, "isVisible", "Enable"),
                new StsBooleanFieldBean(StsSection.class, "drawSurface", "Fill"),
                new StsBooleanFieldBean(StsSection.class, "drawGrid", "Grid"),
                new StsBooleanFieldBean(StsSection.class, "drawEdges", "Edges"),
                new StsBooleanFieldBean(StsSection.class, "drawLines", "Lines"),
                new StsColorComboBoxFieldBean(StsSection.class, "stsColor", "Color", currentModel.getSpectrum("Basic").getStsColors())
	        };
            return displayFields;
		}
		else
		{
            if(faultDisplayFields != null) return faultDisplayFields;
	        faultDisplayFields = new StsFieldBean[]
		    {
                new StsBooleanFieldBean(StsSection.class, "isVisible", "Enable"),
                new StsBooleanFieldBean(StsSection.class, "drawSurface", "Fill"),
                new StsBooleanFieldBean(StsSection.class, "drawGrid", "Grid"),
                new StsBooleanFieldBean(StsSection.class, "drawEdges", "Edges"),
                new StsBooleanFieldBean(StsSection.class, "drawLines", "Lines"),
                new StsButtonListFieldBean(StsSection.class, "displaySideName", "Side", sideStrings, false),
                new StsColorComboBoxFieldBean(StsSection.class, "stsColor", "Color", currentModel.getSpectrum("Basic").getStsColors())
	        };
            return faultDisplayFields;
		}
	}

	public StsFieldBean[] getPropertyFields()
	{
		if (type != StsSection.FAULT)
		{
            if(propertyFields != null) return propertyFields;
            propertyFields = new StsFieldBean[]
		    {
		        new StsFloatFieldBean(StsSection.class, "dipAngle", "Dip Angle")
	        };
            return propertyFields;
		}
		else
		{
            if(faultPropertyFields != null) return faultPropertyFields;
	        faultPropertyFields = new StsFieldBean[]
		    {
                new StsFloatFieldBean(StsSection.class, "leftGap", "Left Gap"),
                new StsFloatFieldBean(StsSection.class, "rightGap", "Right Gap"),
                new StsFloatFieldBean(StsSection.class, "dipAngle", "Dip Angle")
	        };
            return faultPropertyFields;
		}
	}

	public Object[] getChildren()
	{
		return new Object[0];
	}

	public boolean anyDependencies()
	{
		return false;
	}

	public StsObjectPanel getObjectPanel()
	{
		if (type != StsSection.FAULT)
		{
			if (objectPanel == null)
			{
				objectPanel = StsObjectPanel.constructor(this, true);
			}
			return objectPanel;
		}
		else
		{
			if (faultObjectPanel == null)
			{
				faultObjectPanel = StsObjectPanel.constructor(this, true);
			}
			return faultObjectPanel;
		}
	}

	public void treeObjectSelected()
	{
		currentModel.getCreateStsClass("com.Sts.DBTypes.StsSection").selected(this);
	}

	public String getLabel()
	{
		return new String("section: " + getName());
	}

	public String labelString()
	{
		StsLine firstLine = (StsLine) lines.getFirst();
		StsLine lastLine = (StsLine) lines.getLast();
		return new String(getLabel() +
						  " from " + firstLine.lineOnSectionLabel() +
						  " to " + lastLine.lineOnSectionLabel());
	}

	public void logMessage()
	{
		logMessage(labelString());
	}
}
