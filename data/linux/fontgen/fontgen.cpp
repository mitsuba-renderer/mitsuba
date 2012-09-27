#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/fstream.h>
#include <ft2build.h>
#include FT_FREETYPE_H
#include FT_GLYPH_H

using namespace mitsuba;

static const uint32_t KFontGlyphScale = 6;
static const uint32_t KFontGlyphSpacing = 5;

struct Glyph {
	Point2 uv;
	Vector2 uvSize;
	Vector2i pxSize;
	int horizontalAdvance;
	int horizontalBearing;
	int verticalBearing;

	inline Glyph() { }
	inline Glyph(Stream *stream) {
		uv = Point2(
			stream->readSingle(),
			stream->readSingle()
		);
		uvSize = Vector2(
			stream->readSingle(),
			stream->readSingle()
		);
		pxSize = Vector2i(stream);
		horizontalAdvance = stream->readInt();
		horizontalBearing = stream->readInt();
		verticalBearing = stream->readInt();
	}

	void serialize(Stream *stream) {
		stream->writeSingle(uv.x);
		stream->writeSingle(uv.y);
		stream->writeSingle(uvSize.x);
		stream->writeSingle(uvSize.y);
		stream->writeInt(pxSize.x);
		stream->writeInt(pxSize.y);
		stream->writeInt(horizontalBearing);
		stream->writeInt(verticalBearing);
		stream->writeInt(horizontalAdvance);
	}
};

int main(int argc, char **argv) {
	Class::staticInitialization();
	Object::staticInitialization();
	Thread::staticInitialization();
	Logger::staticInitialization();

	if (argc < 5) {
		cout << "Syntax: fontgen <font.ttf> <font.png> <font.desc> <pixel size>" << endl;
		return 0;
	}

	FT_Library library;
	if (FT_Init_FreeType(&library))
		SLog(EError, "Error initializing the FreeType library");

	FT_Face face;
	if (FT_New_Face(library, argv[1], 0, &face))
		SLog(EError, "Unable to load font");

	std::string name = FT_Get_Postscript_Name(face);
	cout << "Font name   : " << name << endl;
	cout << "Glyph count : " << face->num_glyphs << endl;

	/* Set the glyph height */
	int fontSize = atoi(argv[4]);
	if (FT_Set_Pixel_Sizes(face, 0, fontSize))
		SLog(EError, "Error while setting character size");

	int8_t *kerningMatrix = new int8_t[256*256];
	uint32_t maxWidth = 0, maxHeight = 0, maxVertBearing = 0;
	cout << endl;
	memset(kerningMatrix, 0, 256*256);

	cout << "Analyzing.." << endl;
	for (int i=0; i<255; ++i) {
//		int index = FT_Get_Char_Index(face, i);
		int index = i;
		if (FT_Load_Char(face, i, FT_LOAD_RENDER))
			continue;

		uint32_t height = face->glyph->metrics.height >> KFontGlyphScale;
		uint32_t width = face->glyph->metrics.width >> KFontGlyphScale;
		uint32_t vertBearing = face->glyph->metrics.horiBearingY >> KFontGlyphScale;

		for (uint32_t o=0; o<256; o++) {
			FT_Vector kerning;
//			int index2 = FT_Get_Char_Index(face, o);
			int index2 = o;
			FT_Get_Kerning(face, index, index2, FT_KERNING_DEFAULT, &kerning);
			kerningMatrix[i + o*256] = kerning.x >> KFontGlyphScale;
		}

		maxWidth = std::max(maxWidth, width);
		maxHeight = std::max(maxWidth, height);
		if (vertBearing <= fontSize)
			maxVertBearing = std::max(maxVertBearing, vertBearing);
	}
	cout << "Max. width  : " << maxWidth << endl;
	cout << "Max. height : " << maxHeight << endl;
	cout << "Max. vbear  : " << maxVertBearing << endl;

	/* 5 pixels spacing to avoid interpolation artefacts */
	maxWidth += KFontGlyphSpacing;
	maxHeight += KFontGlyphSpacing;

	/* Calculate the total required area in square pixels */
	uint32_t area = 256 * maxWidth * maxHeight;
	
	/* Side length of the square bitmap */
	uint32_t finalSide = 0, side = (uint32_t) sqrtf(area);

	/* Safety margin */
	side += maxWidth > maxHeight ? maxWidth : maxHeight;

	/* Round up to a power of two */
	for (int i=0; i<12 && finalSide < side; i++) {
		finalSide = 1 << i;
	}

	if (finalSide < side) 
		SLog(EError, "The requested font bitmap is bigger than 4096x4096!");

	cout << "Final res.  : " << finalSide << endl;

	/* Number of glyph tiles in each row */
	uint32_t xcount = (uint32_t) std::floor((float) finalSide / (float) maxWidth);

	/* Create an bitmap with alpha */
	ref<Bitmap> bitmap = new Bitmap(Bitmap::ELuminance, Bitmap::EUInt8, Vector2i(finalSide));
	bitmap->clear();
	uint8_t *data = bitmap->getUInt8Data();
	ref<FileStream> fs = new FileStream(argv[3], FileStream::ETruncReadWrite);
	fs->setByteOrder(Stream::ENetworkByteOrder);

	for (uint32_t i=0; i<256; i++) {
		//int index = FT_Get_Char_Index(face, i);
		int index = i;
		if (FT_Load_Char(face, index, FT_LOAD_RENDER)) 
			continue;

		/* Calculate the tile index in the texture */
		uint32_t ytile = i / xcount;
		uint32_t xtile = i % xcount;

		/* Calculate the position in the texture */
		uint32_t xpos = xtile * maxWidth;
		uint32_t ypos = ytile * maxHeight;

		uint8_t *buffer = face->glyph->bitmap.buffer;

		uint32_t width = face->glyph->bitmap.width;
		uint32_t height = face->glyph->bitmap.rows;

		for (uint32_t j=0; j<height; j++) {
			uint8_t *dest = &data[((ypos+j) * finalSide + xpos)];

			for (uint32_t o=0; o<width; o++) {
				/* Copy pixel-by pixel */
				*dest++ = *buffer++;
			}
		}

		Glyph glyph;

		/* Calculate the glyph position on the texture */
		glyph.uv = Point2(
			(float) xpos / (float) finalSide,
			(float) ypos / (float) finalSide
		);

		/* Set the glyph size in pixels */
		glyph.pxSize = Vector2i(width, height);

		/* Set the glyph size on the texture */
		glyph.uvSize = Vector2(
			(float) glyph.pxSize.x / (float) finalSide,
			(float) glyph.pxSize.y / (float) finalSide
		);

		FT_Glyph_Metrics metrics = face->glyph->metrics;

		/* Retrieve further glyph metrics */
//		glyph.horizontalBearing = (int) (metrics.horiBearingX / 64.0f);
		glyph.horizontalBearing = (int) face->glyph->bitmap_left;
		glyph.horizontalAdvance = (int) face->glyph->advance.x / 64;
//		glyph.horizontalAdvance = (int) (metrics.horiAdvance / 64.0f);
		glyph.verticalBearing = (int) (metrics.horiBearingY / 64.0f);
		glyph.serialize(fs);
	}
	fs->write(kerningMatrix, sizeof(int8_t)*256*256);
	fs->close();

	fs = new FileStream(argv[2], FileStream::ETruncReadWrite);
	bitmap->write(Bitmap::EPNG, fs);
	fs->close();

	delete[] kerningMatrix;

	FT_Done_Face(face);
	FT_Done_FreeType(library);
	Class::staticShutdown();
	Object::staticShutdown();
	Thread::staticShutdown();
	return 0;
}
