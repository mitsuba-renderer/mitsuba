#include <mitsuba/render/track.h>

MTS_NAMESPACE_BEGIN

AnimatedTransform::AnimatedTransform(Stream *stream) {
	size_t nTracks = (size_t) stream->readUInt();
	for (size_t i=0; i<nTracks; ++i) {
		AbstractAnimationTrack::EType type = 
			(AbstractAnimationTrack::EType) stream->readUInt();
		AbstractAnimationTrack *track = NULL;
		switch (type) {
			case AbstractAnimationTrack::ETranslationX:
			case AbstractAnimationTrack::ETranslationY:
			case AbstractAnimationTrack::ETranslationZ:
			case AbstractAnimationTrack::EScaleX:
			case AbstractAnimationTrack::EScaleY:
			case AbstractAnimationTrack::EScaleZ:
			case AbstractAnimationTrack::ERotationX:
			case AbstractAnimationTrack::ERotationY:
			case AbstractAnimationTrack::ERotationZ:
				track = new FloatTrack(type, stream);
				break;
			case AbstractAnimationTrack::ETranslationXYZ:
			case AbstractAnimationTrack::EScaleXYZ:
				track = new VectorTrack(type, stream);
				break;
			case AbstractAnimationTrack::ERotationQuat:
				track = new QuatTrack(type, stream);
				break;
			default:
				Log(EError, "Encountered an unknown animation track type (%i)!", type);
		}

		track->incRef();
		m_tracks.push_back(track);
	}
}

void AnimatedTransform::addTrack(AbstractAnimationTrack *track) {
	track->incRef();
	m_tracks.push_back(track);
}

AnimatedTransform::~AnimatedTransform() {
	for (size_t i=0; i<m_tracks.size(); ++i)
		m_tracks[i]->decRef();
}

void AnimatedTransform::serialize(Stream *stream) const {
	stream->writeUInt((uint32_t) m_tracks.size());
	for (size_t i=0; i<m_tracks.size(); ++i)
		m_tracks[i]->serialize(stream);
}

/// Compute the transformation at the specified time value
void AnimatedTransform::eval(Float t, Transform &trafo) const {
	Vector translation(0.0f);
	Vector scale(1.0f);
	Quaternion rotation;

	for (size_t i=0; i<m_tracks.size(); ++i) {
		AbstractAnimationTrack *track = m_tracks[i];
		switch (track->getType()) {
			case AbstractAnimationTrack::ETranslationX:
				translation.x = static_cast<FloatTrack *>(track)->eval(t);
				break;
			case AbstractAnimationTrack::ETranslationY: 
				translation.y = static_cast<FloatTrack *>(track)->eval(t);
				break;
			case AbstractAnimationTrack::ETranslationZ: 
				translation.z = static_cast<FloatTrack *>(track)->eval(t);
				break;
			case AbstractAnimationTrack::ETranslationXYZ: 
				translation = static_cast<VectorTrack *>(track)->eval(t);
				break;
			case AbstractAnimationTrack::EScaleX:
				scale.x = static_cast<FloatTrack *>(track)->eval(t);
				break;
			case AbstractAnimationTrack::EScaleY: 
				scale.y = static_cast<FloatTrack *>(track)->eval(t);
				break;
			case AbstractAnimationTrack::EScaleZ: 
				scale.z = static_cast<FloatTrack *>(track)->eval(t);
				break;
			case AbstractAnimationTrack::EScaleXYZ: 
				scale = static_cast<VectorTrack *>(track)->eval(t);
				break;
			case AbstractAnimationTrack::ERotationQuat: 
				rotation = static_cast<QuatTrack *>(track)->eval(t);
				break;
			default:
				Log(EError, "Encountered an unsupported "
					"animation track type: %i!", track->getType());
		}
	}
	trafo = Transform::translate(translation) * 
		rotation.toTransform() *
		Transform::scale(scale);
}

MTS_IMPLEMENT_CLASS(AbstractAnimationTrack, true, Object)
MTS_IMPLEMENT_CLASS(AnimatedTransform, false, Object)
MTS_NAMESPACE_END
