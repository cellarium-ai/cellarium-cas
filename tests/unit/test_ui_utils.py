from parameterized import parameterized

from cellarium.cas.visualization.ui_utils import ConfigValue


class TestCasClient:
    @parameterized.expand(
        [
            (1),
            (1.1),
            ("a"),
        ]
    )
    def test_config_values_init(self, initValue: any):
        value = ConfigValue(initValue)

        # Make sure simple get works
        assert value.get() == initValue
        assert value.is_dirty() is False
        assert value.get(dirty_read=True) == initValue

    @parameterized.expand(
        [
            (1),
            (1.1),
            ("a"),
        ]
    )
    def test_config_values_simple_rollback(self, initValue: any):
        value = ConfigValue(initValue)

        value.rollback()
        assert value.get() == initValue
        assert value.is_dirty() is False
        assert value.get(dirty_read=True) == initValue

    @parameterized.expand(
        [
            (1, 2),
            (1.1, 2.2),
            ("a", "b"),
        ]
    )
    def test_config_values_set_with_rollback(self, initValue: any, newValue: any):
        value = ConfigValue(initValue)

        # Test pre-commit changes
        value.set(newValue)
        assert value.get() == initValue
        assert value.is_dirty() is True
        assert value.get(dirty_read=True) == newValue

        # Test rollback
        value.rollback()
        assert value.get() == initValue
        assert value.is_dirty() is False
        assert value.get(dirty_read=True) == initValue

    @parameterized.expand(
        [
            (1),
            (1.1),
            ("a"),
        ]
    )
    def test_config_values_simple_commit(self, initValue: any):
        value = ConfigValue(initValue)

        value.commit()
        assert value.get() == initValue
        assert value.is_dirty() is False
        assert value.get(dirty_read=True) == initValue

    @parameterized.expand(
        [
            (1, 2),
            (1.1, 2.2),
            ("a", "b"),
        ]
    )
    def test_config_values_set_with_commit(self, initValue: any, newValue: any):
        value = ConfigValue(initValue)

        # Test pre-commit changes
        value.set(newValue)
        assert value.get() == initValue
        assert value.is_dirty() is True
        assert value.get(dirty_read=True) == newValue

        # Test commit
        value.commit()
        assert value.get() == newValue
        assert value.is_dirty() is False
        assert value.get(dirty_read=True) == newValue

    @parameterized.expand(
        [
            (1, 2),
            (1.1, 2.2),
            ("a", "b"),
        ]
    )
    def test_config_values_reset(self, initValue: any, newValue: any):
        value = ConfigValue(initValue)
        value.set(newValue)
        value.commit()

        value.reset()
        assert value.get() == initValue
        assert value.is_dirty() is False
        assert value.get(dirty_read=True) == initValue

    @parameterized.expand(
        [
            (1, 2, 3),
            (1.1, 2.2, 3.3),
            ("a", "b", "c"),
        ]
    )
    def test_config_values_reset_when_dirty(self, initValue: any, newValue: any, dirtyValue: any):
        value = ConfigValue(initValue)
        value.set(newValue)
        value.commit()
        value.set(dirtyValue)

        value.reset()
        assert value.get() == initValue
        assert value.is_dirty() is False
        assert value.get(dirty_read=True) == initValue
