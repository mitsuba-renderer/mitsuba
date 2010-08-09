// This file is generated, don't edit it!!

#if !defined(XERCESC_INCLUDE_GUARD_ERRHEADER_XMLErrs)
#define XERCESC_INCLUDE_GUARD_ERRHEADER_XMLErrs

#include <xercesc/framework/XMLErrorReporter.hpp>
#include <xercesc/util/XercesDefs.hpp>
#include <xercesc/dom/DOMError.hpp>

XERCES_CPP_NAMESPACE_BEGIN

class XMLErrs
{
public :
    enum Codes
    {
        NoError                            = 0
      , W_LowBounds                        = 1
      , NotationAlreadyExists              = 2
      , AttListAlreadyExists               = 3
      , ContradictoryEncoding              = 4
      , UndeclaredElemInCM                 = 5
      , UndeclaredElemInAttList            = 6
      , XMLException_Warning               = 7
      , XIncludeResourceErrorWarning       = 8
      , XIncludeCannotOpenFile             = 9
      , XIncludeIncludeFailedResourceError   = 10
      , W_HighBounds                       = 11
      , E_LowBounds                        = 12
      , FeatureUnsupported                 = 13
      , TopLevelNoNameComplexType          = 14
      , TopLevelNoNameAttribute            = 15
      , NoNameRefAttribute                 = 16
      , NoNameRefElement                   = 17
      , NoNameRefGroup                     = 18
      , NoNameRefAttGroup                  = 19
      , AnonComplexTypeWithName            = 20
      , AnonSimpleTypeWithName             = 21
      , InvalidElementContent              = 22
      , SimpleTypeContentError             = 23
      , ExpectedSimpleTypeInList           = 24
      , ListUnionRestrictionError          = 25
      , SimpleTypeDerivationByListError    = 26
      , ExpectedSimpleTypeInRestriction    = 27
      , DuplicateFacet                     = 28
      , ExpectedSimpleTypeInUnion          = 29
      , EmptySimpleTypeContent             = 30
      , InvalidSimpleContent               = 31
      , UnspecifiedBase                    = 32
      , InvalidComplexContent              = 33
      , SchemaElementContentError          = 34
      , ContentError                       = 35
      , UnknownSimpleType                  = 36
      , UnknownComplexType                 = 37
      , UnresolvedPrefix                   = 38
      , RefElementNotFound                 = 39
      , TypeNotFound                       = 40
      , TopLevelAttributeNotFound          = 41
      , InvalidChildInComplexType          = 42
      , BaseTypeNotFound                   = 43
      , DatatypeValidatorCreationError     = 44
      , InvalidChildFollowingSimpleContent   = 45
      , InvalidChildFollowingConplexContent   = 46
      , AttributeDefaultFixedValue         = 47
      , NotOptionalDefaultAttValue         = 48
      , DuplicateAttribute                 = 49
      , AttributeWithTypeAndSimpleType     = 50
      , AttributeSimpleTypeNotFound        = 51
      , ElementWithFixedAndDefault         = 52
      , InvalidDeclarationName             = 53
      , ElementWithTypeAndAnonType         = 54
      , NotSimpleOrMixedElement            = 55
      , DisallowedSimpleTypeExtension      = 56
      , InvalidSimpleContentBase           = 57
      , InvalidComplexTypeBase             = 58
      , InvalidChildInSimpleContent        = 59
      , InvalidChildInComplexContent       = 60
      , AnnotationError                    = 61
      , DisallowedBaseDerivation           = 62
      , SubstitutionRepeated               = 63
      , UnionRepeated                      = 64
      , ExtensionRepeated                  = 65
      , ListRepeated                       = 66
      , RestrictionRepeated                = 67
      , InvalidBlockValue                  = 68
      , InvalidFinalValue                  = 69
      , InvalidSubstitutionGroupElement    = 70
      , SubstitutionGroupTypeMismatch      = 71
      , DuplicateElementDeclaration        = 72
      , InvalidAttValue                    = 73
      , AttributeRefContentError           = 74
      , DuplicateRefAttribute              = 75
      , ForbiddenDerivationByRestriction   = 76
      , ForbiddenDerivationByExtension     = 77
      , BaseNotComplexType                 = 78
      , ImportNamespaceDifference          = 79
      , DeclarationNoSchemaLocation        = 80
      , IncludeNamespaceDifference         = 81
      , OnlyAnnotationExpected             = 82
      , InvalidAttributeContent            = 83
      , AttributeRequiredGlobal            = 84
      , AttributeRequiredLocal             = 85
      , AttributeDisallowedGlobal          = 86
      , AttributeDisallowedLocal           = 87
      , InvalidMin2MaxOccurs               = 88
      , AnyAttributeContentError           = 89
      , NoNameGlobalElement                = 90
      , NoCircularDefinition               = 91
      , DuplicateGlobalType                = 92
      , DuplicateGlobalDeclaration         = 93
      , WS_CollapseExpected                = 94
      , Import_1_1                         = 95
      , Import_1_2                         = 96
      , ElemIDValueConstraint              = 97
      , NoNotationType                     = 98
      , EmptiableMixedContent              = 99
      , EmptyComplexRestrictionDerivation   = 100
      , MixedOrElementOnly                 = 101
      , InvalidContentRestriction          = 102
      , ForbiddenDerivation                = 103
      , AtomicItemType                     = 104
      , GroupContentError                  = 105
      , AttGroupContentError               = 106
      , MinMaxOnGroupChild                 = 107
      , DeclarationNotFound                = 108
      , AllContentLimited                  = 109
      , BadMinMaxAllCT                     = 110
      , BadMinMaxAllElem                   = 111
      , DuplicateAttInDerivation           = 112
      , NotExpressibleWildCardIntersection   = 113
      , BadAttDerivation_1                 = 114
      , BadAttDerivation_2                 = 115
      , BadAttDerivation_3                 = 116
      , BadAttDerivation_4                 = 117
      , BadAttDerivation_5                 = 118
      , BadAttDerivation_6                 = 119
      , BadAttDerivation_7                 = 120
      , BadAttDerivation_8                 = 121
      , BadAttDerivation_9                 = 122
      , AllContentError                    = 123
      , RedefineNamespaceDifference        = 124
      , Redefine_InvalidSimpleType         = 125
      , Redefine_InvalidSimpleTypeBase     = 126
      , Redefine_InvalidComplexType        = 127
      , Redefine_InvalidComplexTypeBase    = 128
      , Redefine_InvalidGroupMinMax        = 129
      , Redefine_DeclarationNotFound       = 130
      , Redefine_GroupRefCount             = 131
      , Redefine_AttGroupRefCount          = 132
      , Redefine_InvalidChild              = 133
      , Notation_DeclNotFound              = 134
      , IC_DuplicateDecl                   = 135
      , IC_BadContent                      = 136
      , IC_KeyRefReferNotFound             = 137
      , IC_KeyRefCardinality               = 138
      , IC_XPathExprMissing                = 139
      , AttUseCorrect                      = 140
      , AttDeclPropCorrect3                = 141
      , AttDeclPropCorrect5                = 142
      , AttGrpPropCorrect3                 = 143
      , InvalidTargetNSValue               = 144
      , XMLException_Error                 = 145
      , InvalidRedefine                    = 146
      , InvalidNSReference                 = 147
      , NotAllContent                      = 148
      , InvalidAnnotationContent           = 149
      , InvalidFacetName                   = 150
      , InvalidXMLSchemaRoot               = 151
      , CircularSubsGroup                  = 152
      , ELTSchemaNS                        = 153
      , InvalidAttTNS                      = 154
      , NSDeclInvalid                      = 155
      , DOMLevel1Node                      = 156
      , E_HighBounds                       = 157
      , F_LowBounds                        = 158
      , EntityExpansionLimitExceeded       = 159
      , ExpectedCommentOrCDATA             = 160
      , ExpectedAttrName                   = 161
      , ExpectedNotationName               = 162
      , NoRepInMixed                       = 163
      , ExpectedDefAttrDecl                = 164
      , ExpectedEqSign                     = 165
      , ExpectedElementName                = 166
      , CommentsMustStartWith              = 167
      , InvalidDocumentStructure           = 168
      , ExpectedDeclString                 = 169
      , BadXMLVersion                      = 170
      , UnsupportedXMLVersion              = 171
      , UnterminatedXMLDecl                = 172
      , BadXMLEncoding                     = 173
      , BadStandalone                      = 174
      , UnterminatedComment                = 175
      , PINameExpected                     = 176
      , UnterminatedPI                     = 177
      , InvalidCharacter                   = 178
      , UnterminatedStartTag               = 179
      , ExpectedAttrValue                  = 180
      , UnterminatedEndTag                 = 181
      , ExpectedAttributeType              = 182
      , ExpectedEndOfTagX                  = 183
      , ExpectedMarkup                     = 184
      , NotValidAfterContent               = 185
      , ExpectedComment                    = 186
      , ExpectedCommentOrPI                = 187
      , ExpectedWhitespace                 = 188
      , NoRootElemInDOCTYPE                = 189
      , ExpectedQuotedString               = 190
      , ExpectedPublicId                   = 191
      , InvalidPublicIdChar                = 192
      , UnterminatedDOCTYPE                = 193
      , InvalidCharacterInIntSubset        = 194
      , UnexpectedWhitespace               = 195
      , InvalidCharacterInAttrValue        = 196
      , ExpectedMarkupDecl                 = 197
      , TextDeclNotLegalHere               = 198
      , ConditionalSectInIntSubset         = 199
      , ExpectedPEName                     = 200
      , UnterminatedEntityDecl             = 201
      , InvalidCharacterRef                = 202
      , UnterminatedCharRef                = 203
      , ExpectedEntityRefName              = 204
      , EntityNotFound                     = 205
      , NoUnparsedEntityRefs               = 206
      , UnterminatedEntityRef              = 207
      , RecursiveEntity                    = 208
      , PartialMarkupInEntity              = 209
      , UnterminatedElementDecl            = 210
      , ExpectedContentSpecExpr            = 211
      , ExpectedAsterisk                   = 212
      , UnterminatedContentModel           = 213
      , ExpectedSystemOrPublicId           = 214
      , UnterminatedNotationDecl           = 215
      , ExpectedSeqChoiceLeaf              = 216
      , ExpectedChoiceOrCloseParen         = 217
      , ExpectedSeqOrCloseParen            = 218
      , ExpectedEnumValue                  = 219
      , ExpectedEnumSepOrParen             = 220
      , UnterminatedEntityLiteral          = 221
      , MoreEndThanStartTags               = 222
      , ExpectedOpenParen                  = 223
      , AttrAlreadyUsedInSTag              = 224
      , BracketInAttrValue                 = 225
      , Expected2ndSurrogateChar           = 226
      , ExpectedEndOfConditional           = 227
      , ExpectedIncOrIgn                   = 228
      , ExpectedINCLUDEBracket             = 229
      , UnexpectedEOE                      = 230
      , PEPropogated                       = 231
      , ExtraCloseSquare                   = 232
      , PERefInMarkupInIntSubset           = 233
      , EntityPropogated                   = 234
      , ExpectedNumericalCharRef           = 235
      , ExpectedOpenSquareBracket          = 236
      , BadSequenceInCharData              = 237
      , IllegalSequenceInComment           = 238
      , UnterminatedCDATASection           = 239
      , ExpectedNDATA                      = 240
      , NDATANotValidForPE                 = 241
      , HexRadixMustBeLowerCase            = 242
      , DeclStringRep                      = 243
      , DeclStringsInWrongOrder            = 244
      , NoExtRefsInAttValue                = 245
      , XMLDeclMustBeLowerCase             = 246
      , ExpectedEntityValue                = 247
      , BadDigitForRadix                   = 248
      , EndedWithTagsOnStack               = 249
      , NestedCDATA                        = 250
      , UnknownPrefix                      = 251
      , PartialTagMarkupError              = 252
      , EmptyMainEntity                    = 253
      , CDATAOutsideOfContent              = 254
      , Unexpected2ndSurrogateChar         = 255
      , NoPIStartsWithXML                  = 256
      , XMLDeclMustBeFirst                 = 257
      , XMLVersionRequired                 = 258
      , StandaloneNotLegal                 = 259
      , EncodingRequired                   = 260
      , ColonNotLegalWithNS                = 261
      , XMLException_Fatal                 = 262
      , BadSchemaLocation                  = 263
      , SchemaScanFatalError               = 264
      , IllegalRefInStandalone             = 265
      , PEBetweenDecl                      = 266
      , NoEmptyStrNamespace                = 267
      , NoUseOfxmlnsAsPrefix               = 268
      , NoUseOfxmlnsURI                    = 269
      , PrefixXMLNotMatchXMLURI            = 270
      , XMLURINotMatchXMLPrefix            = 271
      , NoXMLNSAsElementPrefix             = 272
      , CT_SimpleTypeChildRequired         = 273
      , InvalidRootElemInDOCTYPE           = 274
      , InvalidElementName                 = 275
      , InvalidAttrName                    = 276
      , InvalidEntityRefName               = 277
      , DuplicateDocTypeDecl               = 278
      , XIncludeOrphanFallback             = 279
      , XIncludeNoHref                     = 280
      , XIncludeXPointerNotSupported       = 281
      , XIncludeInvalidParseVal            = 282
      , XIncludeMultipleFallbackElems      = 283
      , XIncludeIncludeFailedNoFallback    = 284
      , XIncludeCircularInclusionLoop      = 285
      , XIncludeCircularInclusionDocIncludesSelf   = 286
      , XIncludeDisallowedChild            = 287
      , XIncludeConflictingNotation        = 288
      , XIncludeConflictingEntity          = 289
      , F_HighBounds                       = 290
    };

    static bool isFatal(const XMLErrs::Codes toCheck)
    {
        return ((toCheck >= F_LowBounds) && (toCheck <= F_HighBounds));
    }

    static bool isWarning(const XMLErrs::Codes toCheck)
    {
        return ((toCheck >= W_LowBounds) && (toCheck <= W_HighBounds));
    }

    static bool isError(const XMLErrs::Codes toCheck)
    {
        return ((toCheck >= E_LowBounds) && (toCheck <= E_HighBounds));
    }

    static XMLErrorReporter::ErrTypes errorType(const XMLErrs::Codes toCheck)
    {
       if ((toCheck >= W_LowBounds) && (toCheck <= W_HighBounds))
           return XMLErrorReporter::ErrType_Warning;
       else if ((toCheck >= F_LowBounds) && (toCheck <= F_HighBounds))
            return XMLErrorReporter::ErrType_Fatal;
       else if ((toCheck >= E_LowBounds) && (toCheck <= E_HighBounds))
            return XMLErrorReporter::ErrType_Error;
       return XMLErrorReporter::ErrTypes_Unknown;
    }
    static DOMError::ErrorSeverity  DOMErrorType(const XMLErrs::Codes toCheck)
    {
       if ((toCheck >= W_LowBounds) && (toCheck <= W_HighBounds))
           return DOMError::DOM_SEVERITY_WARNING;
       else if ((toCheck >= F_LowBounds) && (toCheck <= F_HighBounds))
            return DOMError::DOM_SEVERITY_FATAL_ERROR;
       else return DOMError::DOM_SEVERITY_ERROR;
    }

private:
    // -----------------------------------------------------------------------
    //  Unimplemented constructors and operators
    // -----------------------------------------------------------------------
    XMLErrs();
};

XERCES_CPP_NAMESPACE_END

#endif

