#include <mitsuba/core/track.h>
#include <mitsuba/core/aabb.h>
#include <Eigen/SVD>
#include <Eigen/Geometry>

MTS_NAMESPACE_BEGIN

AnimatedTransform::AnimatedTransform(const AnimatedTransform *trafo)
		: m_transform(trafo->m_transform) {
	m_tracks.reserve(trafo->getTrackCount());
	for (size_t i=0; i<trafo->getTrackCount(); ++i) {
		AbstractAnimationTrack *track = trafo->getTrack(i)->clone();
		m_tracks.push_back(track);
		track->incRef();
	}
}

AnimatedTransform::AnimatedTransform(Stream *stream) {
	size_t nTracks = stream->readSize();
	if (nTracks == 0) {
		m_transform = Transform(stream);
	} else {
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
}

void AnimatedTransform::addTrack(AbstractAnimationTrack *track) {
	track->incRef();
	m_tracks.push_back(track);
}

AABB1 AnimatedTransform::getTimeBounds() const {
	if (m_tracks.size() == 0)
		return AABB1(0.0f, 0.0f);

	Float min =  std::numeric_limits<Float>::infinity();
	Float max = -std::numeric_limits<Float>::infinity();

	for (size_t i=0; i<m_tracks.size(); ++i) {
		const AbstractAnimationTrack *track = m_tracks[i];
		size_t size = track->getSize();
		SAssert(size > 0);
		min = std::min(min, track->getTime(0));
		max = std::max(max, track->getTime(size-1));
	}

	return AABB1(min, max);
}

AABB AnimatedTransform::getTranslationBounds() const {
	if (m_tracks.size() == 0) {
		Point p = m_transform(Point(0.0f));
		return AABB(p, p);
	}

	AABB aabb;

	for (size_t i=0; i<m_tracks.size(); ++i) {
		const AbstractAnimationTrack *absTrack = m_tracks[i];
		switch (absTrack->getType()) {
			case AbstractAnimationTrack::ETranslationX:
			case AbstractAnimationTrack::ETranslationY:
			case AbstractAnimationTrack::ETranslationZ: {
					int idx  = absTrack->getType() - AbstractAnimationTrack::ETranslationX;
					const FloatTrack *track =
						static_cast<const FloatTrack *>(absTrack);
					for (size_t j=0; j<track->getSize(); ++j) {
						Float value = track->getValue(j);
						aabb.max[idx] = std::max(aabb.max[idx], value);
						aabb.min[idx] = std::min(aabb.min[idx], value);
					}
				}
				break;

			case AbstractAnimationTrack::ETranslationXYZ: {
					const VectorTrack *track =
						static_cast<const VectorTrack *>(absTrack);
					for (size_t j=0; j<track->getSize(); ++j)
						aabb.expandBy(Point(track->getValue(j)));
				}
				break;
			default:
				break;
		}
	}
	for (int i=0; i<3; ++i) {
		if (aabb.min[i] > aabb.max[i])
			aabb.min[i] = aabb.max[i] = 0.0f;
	}

	return aabb;
}

AABB AnimatedTransform::getSpatialBounds(const AABB &aabb) const {
	AABB result;

	if (m_tracks.size() == 0) {
		for (int j=0; j<8; ++j)
			result.expandBy(m_transform(aabb.getCorner(j)));
	} else {
		/* Compute approximate bounds */
		int nSteps = 100;
		AABB1 timeBounds = getTimeBounds();
		Float step = timeBounds.getExtents().x / (nSteps-1);

		for (int i=0; i<nSteps; ++i) {
			const Transform &trafo = eval(timeBounds.min.x + step * i);
			for (int j=0; j<8; ++j)
				result.expandBy(trafo(aabb.getCorner(j)));
		}
	}

	return result;
}

AnimatedTransform::~AnimatedTransform() {
	for (size_t i=0; i<m_tracks.size(); ++i)
		m_tracks[i]->decRef();
}

void AnimatedTransform::sortAndSimplify() {
	bool isStatic = true;

	for (size_t i=0; i<m_tracks.size(); ++i) {
		AbstractAnimationTrack *track = m_tracks[i];
		bool isNeeded = false;
		switch (track->getType()) {
			case AbstractAnimationTrack::ETranslationX:
			case AbstractAnimationTrack::ETranslationY:
			case AbstractAnimationTrack::ETranslationZ:
			case AbstractAnimationTrack::ERotationX:
			case AbstractAnimationTrack::ERotationY:
			case AbstractAnimationTrack::ERotationZ:
			case AbstractAnimationTrack::EScaleX:
			case AbstractAnimationTrack::EScaleY:
			case AbstractAnimationTrack::EScaleZ:
				isNeeded = static_cast<FloatTrack *>(track)->sortAndSimplify();
				break;
			case AbstractAnimationTrack::ETranslationXYZ:
			case AbstractAnimationTrack::EScaleXYZ:
				isNeeded = static_cast<VectorTrack *>(track)->sortAndSimplify();
				break;
			case AbstractAnimationTrack::ERotationQuat:
				isNeeded = static_cast<QuatTrack *>(track)->sortAndSimplify();
				break;
			default:
				Log(EError, "Encountered an unsupported "
					"animation track type: %i!", track->getType());
		}
		if (isNeeded) {
			isStatic &= track->getSize() == 1;
		} else {
			m_tracks.erase(m_tracks.begin() + i);
			track->decRef();
			--i;
		}
	}

	if (isStatic) {
		Transform temp;
		temp = eval(0);
		m_transform = temp;
		for (size_t i=0; i<m_tracks.size(); ++i)
			m_tracks[i]->decRef();
		m_tracks.clear();
	}
}


const AbstractAnimationTrack *AnimatedTransform::findTrack(AbstractAnimationTrack::EType type) const {
	for (size_t i=0; i<m_tracks.size(); ++i) {
		AbstractAnimationTrack *track = m_tracks[i];
		if (track->getType() == type)
			return track;
	}
	return NULL;
}
AbstractAnimationTrack *AnimatedTransform::findTrack(AbstractAnimationTrack::EType type) {
	for (size_t i=0; i<m_tracks.size(); ++i) {
		AbstractAnimationTrack *track = m_tracks[i];
		if (track->getType() == type)
			return track;
	}
	return NULL;
}

void AnimatedTransform::prependScale(const Vector &scale) {
	FloatTrack *trackX = (FloatTrack *) findTrack(AbstractAnimationTrack::EScaleX);
	FloatTrack *trackY = (FloatTrack *) findTrack(AbstractAnimationTrack::EScaleY);
	FloatTrack *trackZ = (FloatTrack *) findTrack(AbstractAnimationTrack::EScaleZ);
	VectorTrack *trackXYZ = (VectorTrack *) findTrack(AbstractAnimationTrack::EScaleXYZ);

	if (m_tracks.empty()) {
		m_transform = m_transform * Transform::scale(scale);
	} else if (trackXYZ) {
		trackXYZ->prependTransformation(scale);
	} else if (trackX && trackY && trackZ) {
		if (trackX) {
			trackX->prependTransformation(scale.x);
		} else {
			trackX = new FloatTrack(AbstractAnimationTrack::EScaleX);
			trackX->append(0.0f, scale.x); addTrack(trackX);
		}

		if (trackY) {
			trackY->prependTransformation(scale.y);
		} else {
			trackY = new FloatTrack(AbstractAnimationTrack::EScaleY);
			trackY->append(0.0f, scale.y); addTrack(trackY);
		}

		if (trackZ) {
			trackZ->prependTransformation(scale.z);
		} else {
			trackZ = new FloatTrack(AbstractAnimationTrack::EScaleZ);
			trackZ->append(0.0f, scale.z); addTrack(trackZ);
		}
	} else {
		trackXYZ = new VectorTrack(AbstractAnimationTrack::EScaleXYZ);
		trackXYZ->append(0.0f, scale);
		addTrack(trackXYZ);
	}
}

void AnimatedTransform::collectKeyframes(std::set<Float> &result) const {
	for (size_t i=0; i<m_tracks.size(); ++i) {
		const AbstractAnimationTrack *track = m_tracks[i];

		for (size_t j=0; j<track->getSize(); ++j)
			result.insert(track->getTime(j));
	}

	if (result.size() == 0)
		result.insert((Float) 0);
}

void AnimatedTransform::serialize(Stream *stream) const {
	stream->writeSize(m_tracks.size());
	if (m_tracks.size() == 0) {
		m_transform.serialize(stream);
	} else {
		for (size_t i=0; i<m_tracks.size(); ++i)
			m_tracks[i]->serialize(stream);
	}
}

void AnimatedTransform::TransformFunctor::operator()(const Float &t, Transform &trafo) const {
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

	trafo = Transform::translate(translation);

	if (!rotation.isIdentity())
		trafo = trafo * rotation.toTransform();

	if (scale != Vector(0.0f))
		trafo = trafo * Transform::scale(scale);
}

void AnimatedTransform::appendTransform(Float time, const Transform &trafo) {
	/* Compute the polar decomposition and insert into the animated transform;
	   uh oh.. we have to get rid of the two separate matrix libraries at some point :) */
	typedef Eigen::Matrix<Float, 3, 3> EMatrix;

	if (m_tracks.size() == 0) {
		ref<VectorTrack> translation = new VectorTrack(VectorTrack::ETranslationXYZ);
		ref<QuatTrack> rotation = new QuatTrack(VectorTrack::ERotationQuat);
		ref<VectorTrack> scaling = new VectorTrack(VectorTrack::EScaleXYZ);
		translation->reserve(2);
		rotation->reserve(2);
		scaling->reserve(2);
		addTrack(translation);
		addTrack(rotation);
		addTrack(scaling);
	} else if (m_tracks.size() != 3 ||
			m_tracks[0]->getType() != VectorTrack::ETranslationXYZ ||
			m_tracks[1]->getType() != VectorTrack::ERotationQuat ||
			m_tracks[2]->getType() != VectorTrack::EScaleXYZ) {
		Log(EError, "AnimatedTransform::appendTransform(): unsupported internal configuration!");
	}

	const Matrix4x4 m = trafo.getMatrix();
	EMatrix A;

	A << m(0, 0), m(0, 1), m(0, 2),
		 m(1, 0), m(1, 1), m(1, 2),
		 m(2, 0), m(2, 1), m(2, 2);

	Eigen::JacobiSVD<EMatrix> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
	EMatrix U = svd.matrixU(), V = svd.matrixV(), S = svd.singularValues().asDiagonal();

	EMatrix Q = U*V.transpose();
	EMatrix P = V*S*V.transpose();

	if (Q.determinant() < 0) {
		Q = -Q; P = -P;
	}

	VectorTrack *translation = (VectorTrack *) m_tracks[0];
	QuatTrack *rotation = (QuatTrack *) m_tracks[1];
	VectorTrack *scaling = (VectorTrack *) m_tracks[2];

	rotation->append(time, Quaternion::fromMatrix(
		Matrix4x4(
			Q(0, 0), Q(0, 1), Q(0, 2), 0.0f,
			Q(1, 0), Q(1, 1), Q(1, 2), 0.0f,
			Q(2, 0), Q(2, 1), Q(2, 2), 0.0f,
			0.0f,    0.0f,    0.0f,    1.0f
		)
	));

	scaling->append(time, Vector(P(0, 0), P(1, 1), P(2, 2)));
	translation->append(time, Vector(m(0, 3), m(1, 3), m(2, 3)));
}

std::string AnimatedTransform::toString() const {
	if (m_tracks.size() == 0) {
		return m_transform.toString();
	} else {
		std::ostringstream oss;
		oss << "AnimatedTransform[tracks=" << m_tracks.size() << "]";
		return oss.str();
	}
}

MTS_IMPLEMENT_CLASS(AbstractAnimationTrack, true, Object)
MTS_IMPLEMENT_CLASS(AnimatedTransform, false, Object)
MTS_NAMESPACE_END
